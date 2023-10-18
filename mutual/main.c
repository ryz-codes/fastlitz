/* main.c
 * Matlab gateway file with inputs as filaments
 * Defined according to O, L, W, H
 *************************************************************/

#include "mutual.h"

/*---------------------------------------------------------------------------------------- */
/* Prototypes */
/*---------------------------------------------------------------------------------------- */
void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]);
int convFils(FILAMENT* fil1, FILAMENT* fil2, double *O, double *L, double *W, double *H);
#define CONVFILS_SELF TRUE
#define CONVFILS_MUTUAL FALSE

/*---------------------------------------------------------------------------------------- */
/* MATLAB MEX INTERFACE */
/*---------------------------------------------------------------------------------------- */
#define XX0 0
#define YY0 1
#define ZZ0 2
#define XX1 3
#define YY1 4
#define ZZ1 5

void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{   
	double *O, *L, *W, *H;
    double *Oin, *Lin, *Win, *Hin;
    double *out; /*, *stats;*/
    int ind;
    /*int dispstats = TRUE;*/
    size_t n;
    FILAMENT fil1, fil2;
    
    /*if (nlhs==2)
        dispstats = FALSE;*/
    if (nlhs>2)
        mexErrMsgTxt("Wrong number of output parameters, usage:  M = mutual(O, L, W, H)");
    if (nrhs!=4)
        mexErrMsgTxt("Wrong number of input parameters, usage:  M = mutual(O, L, W, H)");
    if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]) || !mxIsDouble(prhs[3]))
        mexErrMsgTxt("mutual: Input arguments must be double.");

    /* Prepare input */
    Oin = mxGetPr(prhs[0]); /* O - origin */
    Lin = mxGetPr(prhs[1]); /* L - length */
    Win = mxGetPr(prhs[2]); /* W - width */
    Hin = mxGetPr(prhs[3]); /* H - height */
    
    /* number of elements */
    n = mxGetN(prhs[0]);
    
    /* Error check for m == 6 */
    if(!(6 == mxGetM(prhs[0]) &&
    		 6 == mxGetM(prhs[1]) &&
    		 6 == mxGetM(prhs[2]) &&
    		 6 == mxGetM(prhs[3])))
    		 mexErrMsgTxt("mutual: Must have six rows for each of the four inputs");
    if(n != mxGetN(prhs[1]) ||
       n != mxGetN(prhs[2]) ||
       n != mxGetN(prhs[3]))
    		mexErrMsgTxt("mutual: Must have same number of columns");   		   
    
    /* Output */
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    out = mxGetPr(plhs[0]);
    
    /*if (!dispstats) {
        plhs[1] = mxCreateDoubleMatrix(6, 1, mxREAL);
        stats = mxGetPr(plhs[1]);
    }*/
        
            
    /* Reset counters 
    num_exact_mutual=0;
    num_fourfil=0;
    num_mutualfil=0;
    num_quadFil=0;
    num_perp=0;
    num_self=0;*/
    
    #ifdef _OPENMP
        omp_set_num_threads(omp_get_num_procs());  
    #endif
    
    #pragma omp parallel for default(none) shared(out,n,Oin,Lin,Win,Hin) private(O,L,W,H,fil1,fil2,ind)
    for (ind=0;ind<n;ind++) {
          /* Increment pointers */
          O = Oin + 6*ind; L = Lin + 6*ind; W = Win + 6*ind; H = Hin + 6*ind;
          
          /* Translate and Output */
          if (convFils(&fil1,&fil2,O,L,W,H)) {
              out[ind] = self(fil1.width, fil1.length, fil1.height); 
          }
          else {
              out[ind] = mutual(&fil1,&fil2);   
          }
            
    }
    /*More debug tools */ 
    /*
    if (dispstats) {
         mexPrintf("Calls to exact_mutual: %15d\n",num_exact_mutual); 
         mexPrintf("         cuber:        %15d\n",num_quadFil); 
         mexPrintf("         fourfils:     %15d\n",num_fourfil); 
         mexPrintf("         mutualfil:    %15d\n",num_mutualfil); 
         mexPrintf("Number perpendicular:  %15d\n",num_perp); 
         mexPrintf("Number self terms:     %15d\n",num_self); 
         mexPrintf("\n"); 
    } else {
        stats[0] = num_exact_mutual;
        stats[1] = num_quadFil;
        stats[2] = num_fourfil;
        stats[3] = num_mutualfil;
        stats[4] = num_perp;
        stats[5] = num_self;
    }
    return;
     */
}

/*----------------------------------------------------------------------------------------
 * Interface to Matt's filament definition
 *---------------------------------------------------------------------------------------- */

int convFils(FILAMENT* fil1, FILAMENT* fil2, double *O, double *L, double *W, double *H) {
      double P[6], Q[6];
      int ii;
      
      int chkslf = (1 == 1);
      for (ii=0; ii<3; ii++) {
          chkslf = chkslf && nearzero(O[ii]-O[ii+3]);
          chkslf = chkslf && nearzero(L[ii]-L[ii+3]);
          chkslf = chkslf && nearzero(W[ii]-W[ii+3]);
          chkslf = chkslf && nearzero(H[ii]-H[ii+3]);
      }

      /* Convert to the FILAMENT format         */
      /* Do fil1 first. */
      /* length, width, height */
      fil1->length = mag(L[XX0],L[YY0],L[ZZ0]);
      fil1->width = mag(W[XX0],W[YY0],W[ZZ0]);
      fil1->height = mag(H[XX0],H[YY0],H[ZZ0]);

      if(chkslf) {
         return CONVFILS_SELF;
      } else {
          /* length, width, height */
          fil2->length = mag(L[XX1],L[YY1],L[ZZ1]);
          fil2->width = mag(W[XX1],W[YY1],W[ZZ1]);
          fil2->height = mag(H[XX1],H[YY1],H[ZZ1]);

          /* Make end points.  */
          /* Note that FastHenry defines the origin at the center of the  */
          /* square, but the OLWH notation defines the origin at the corner. */
          /* We define P to be the FastHenry origin and Q to be the end point. */
          for (ii = 0; ii < 6; ii++) {
              P[ii] = O[ii] + W[ii]/2 + H[ii]/2;  
              Q[ii] = O[ii] + L[ii] + W[ii]/2 + H[ii]/2;
          }
          fil1->x[0] = P[XX0]; fil1->y[0] = P[YY0]; fil1->z[0] = P[ZZ0];
          fil1->x[1] = Q[XX0]; fil1->y[1] = Q[YY0]; fil1->z[1] = Q[ZZ0];
          fil2->x[0] = P[XX1]; fil2->y[0] = P[YY1]; fil2->z[0] = P[ZZ1];
          fil2->x[1] = Q[XX1]; fil2->y[1] = Q[YY1]; fil2->z[1] = Q[ZZ1];

          /* lenvect */
          fil1->lenvect = L; fil2->lenvect = L+3;

          /* widvect */
          fil1->widvect = W; fil2->widvect = W+3;
          
          return  CONVFILS_MUTUAL;
      }
    
}