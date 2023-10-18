/* CUBER.C
 * Implements the freespace inductive 1/R volume to volume interactions as a 
 * collection of surface-surface interactions, each governed by a Green's 
 * function of -R/2. 
 * 
 * The dimensions of the source/target boxes are respective labeled:
 *    u = [len1 wid1 hei1], v = [len2 wid2 hei2].
 * The bases of the source/target boxes are labeled:
 *    U = [ux uy uz], V = [vx vy vz].
 *
 * surf_integral() evaluates a single surface to surface interaction, with 
 * one surface lying flat on the xy plane. This uses the sparse grid 
 * quadrature nodes stored in KPU_2_3.h
 *
 * sixface() evaluates the interaction between one surface with six target 
 * surfaces, defined by the transformed basis matrix U' * V = TV, the 
 * dimensions vector v, the bottom-left corner and the top-right corner.
 *
 * cuber() performs the actual calculation of volume-to-volume interactions
 * by rotating everything to the U basis.
 *
 * Version: 1.0.1 (April 14th, 2013)
 * Author: Richard Y. Zhang (ryz@mit.edu)
 *
 */

#include "mutual.h"
#include "KPU_2_2.h"

#define ONE3RD 0.333333333333333333
#define ONE6TH 0.166666666666666666
double surf_integral (double L, double W, double O[3], double L2[3], 
        double W2[3], double len2, double wid2)
{
    double xp[4], yp[4]; /* Definite integration limits */
    double x,y,z,x2,y2,r,z2,z3d3,z2x3,ret[4];
    double tot = 0; 
    int k,i;
    
    /* Quadrature loop */
    for(k=0;k<Q_NUM_WEIGHTS;k++) {
        /* Quadrature points */
        x = O[0] + len2*L2[0]*Qnodes[0][k] + wid2*W2[0]*Qnodes[1][k];
        y = O[1] + len2*L2[1]*Qnodes[0][k] + wid2*W2[1]*Qnodes[1][k];
        z = O[2] + len2*L2[2]*Qnodes[0][k] + wid2*W2[2]*Qnodes[1][k];
        
        /* Set up z-based variables that don't change between def. int*/
        z2 = z*z;
        z3d3 = ONE3RD*z*z2; 
        z2x3 = z2*3; 
        
        /* The definite integral is performed over four indefinite integrals
         * Let us set up the limits ahead of time. 
         *
         * The fact that this is embedded within this function is to avoid 
         * many function calls and the associated overhead */
        xp[0] = x; xp[1] = x-L; xp[2] = xp[0]; xp[3] = xp[1];
        yp[0] = y; yp[1] = yp[0]; yp[2] = y-W; yp[3] = yp[2];

        for(i=0;i<4;i++) {
            x2 = xp[i]*xp[i];
            y2 = yp[i]*yp[i];
            r = sqrt(x2+y2+z2);
            ret[i] = ONE3RD*xp[i]*yp[i]*r;

            /* If clauses to avoid divide by zero*/
            if (y2 > EPS || z2 > EPS)
                ret[i] += ONE6TH * (y2 + z2x3) * yp[i] * log((xp[i]+r)/sqrt(y2+z2));
            if (x2 > EPS || z2 > EPS)
                ret[i] += ONE6TH * (x2 + z2x3) * xp[i] * log((yp[i]+r)/sqrt(x2+z2));

            /* Atan term */   
            if ((z != 0))
                ret[i] -= z3d3 * (atan((xp[i]*yp[i])/(r*z))); 
        }
        /* Decide to add or subtract 
         * Add the first and last, subtract the middle two. */
        tot += (ret[0] - ret[1] - ret[2] + ret[3]) * Qweights[k];
    }
    
    return tot * len2 * wid2;
}

/* Extract the governing directional vectors from the TV matrix
 * We use #define here to avoid unnecessarily creating / destroying pointers
 */
#define L2V (TV[0]+ax)
#define W2V (TV[1]+ax)
#define H2V (TV[2]+ax)
double sixface(int ax, double TV[3][5], double *Db, double *Dt, double len1, 
        double wid1, double *v) {
    
    /* Calculate dot products of the normal of this current plane with the 
     * plane normals of cuboid 2 */
    double dotx = TV[0][ax+2];
    double doty = TV[1][ax+2];
    double dotz = TV[2][ax+2];
    
    double M = 0;
    /* Compute the contributions of each plane if the dot product isnt zero
     */
    if (!nearzero(dotz)) {
        M += dotz * (surf_integral(len1,wid1,Dt,L2V,W2V,-v[0],-v[1]) 
                   - surf_integral(len1,wid1,Db,L2V,W2V,v[0],v[1]));
    }
    if (!nearzero(doty)) {
        M += doty * (surf_integral(len1,wid1,Dt,H2V,L2V,-v[2],-v[0]) 
                   - surf_integral(len1,wid1,Db,H2V,L2V,v[2],v[0]));
    }
    if (!nearzero(dotx)) {
        M += dotx * (surf_integral(len1,wid1,Dt,W2V,H2V,-v[1],-v[2]) 
                   - surf_integral(len1,wid1,Db,W2V,H2V,v[1],v[2]));    
    }
    return M;
}

double cuber(FILAMENT *fil1, FILAMENT *fil2) {
    double u[5], v[5];
    double TV[3][5], U[3][3], V[3][3], O[3];
    double D2b[5], D2t[5];
    double D1b[5] = {0,0,0,0,0};
    double D1t[5] = {0,0,0,0,0};
    
    /* Temporary variables */
    double wid12, wid22, hei12, hei22;
    int i,j; /* column, row*/
    double M = 0;
    
    /* Source dimensions vector*/
    u[0] = fil1->length; u[1] = fil1->width; u[2] = fil1->height; 
    u[3]=u[0]; u[4]=u[1];
    
    /* Source basis matrix*/
    U[0][0] = fil1->lenvect[0] / u[0];
    U[0][1] = fil1->lenvect[1] / u[0];
    U[0][2] = fil1->lenvect[2] / u[0];
    get_wid(fil1,U[1]);
    get_height(fil1,U[1],U[2]);
    
    /* Target dimensions vector*/
    v[0] = fil2->length; v[1] = fil2->width; v[2] = fil2->height; 
    v[3]=v[0]; v[4]=v[1];
    
    /* Target basis matrix*/
    V[0][0] = fil2->lenvect[0] / v[0];
    V[0][1] = fil2->lenvect[1] / v[0];
    V[0][2] = fil2->lenvect[2] / v[0];
    get_wid(fil2,V[1]);
    get_height(fil2,V[1],V[2]);
    
    /* Distance between the origins (center of filaments, not corners) */
    O[0] = fil2->x[0] - fil1->x[0]; 
    O[1] = fil2->y[0] - fil1->y[0]; 
    O[2] = fil2->z[0] - fil1->z[0]; 
    
    /*
     * Set up rotation matrices. See the accompanying matlab files for details
     */
    wid12 = u[1]/2; wid22 = v[1]/2; hei12 = u[2]/2; hei22 = v[2]/2; 
    for (i=0;i<3;i++) { /* Cols of U, V, TV, gives ux, uy, uz etc */
        /* Adjust for offset. After this, O[i] is now dist btwn corners*/
        O[i] += U[1][i]*wid12 - V[1][i]*wid22 + U[2][i]*hei12 - V[2][i]*hei22;
                
        for (j=0;j<3;j++) { /* Rows of U, V, TV, gives inner dimensions */
            /* Form TV*/
            TV[i][j] = U[j][0]*V[i][0] + U[j][1]*V[i][1] + U[j][2]*V[i][2];
            
            /* Form D1b */
            D1b[j] += U[j][i]*O[i];
        }
        
        /* Replicate the rows */
        TV[i][3] = TV[i][0];
        TV[i][4] = TV[i][1];    
        
        /*mexPrintf("TV[%d]: %g %g %g %g %g\n", i, TV[i][0], TV[i][1], TV[i][2], TV[i][3], TV[i][4]);*/
    }
    /* Replicate the rows */
    D1b[3] = D1b[0]; 
    D1b[4] = D1b[1]; 
    
    for (i=0;i<5;i++) {
        /* Form D1t */
        D1t[i] = D1b[i]; 
        for (j=0;j<3;j++)
            D1t[i] += TV[j][i] * v[j];
        
        /* Form D2b and D2t */
        D2b[i] = D1b[i] - u[i];
        D2t[i] = D1t[i] - u[i];
    }
    
    /* Do the faces for each governing axis*/
    for (i=0;i<3;i++) {
        /* Front face*/
        M += sixface(i,TV,D1b+i,D1t+i,u[i],u[i+1],v);
        /* Rear face*/
        M -= sixface(i,TV,D2b+i,D2t+i,-u[i],-u[i+1],v);
    }
    
    /* NB: TV[0][0] = dot(ux,vx) */
    return M * TV[0][0] / (2.0*u[1]*u[2]*v[1]*v[2]);
}