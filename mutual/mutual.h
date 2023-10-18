#ifndef MUTUAL_H
#define MUTUAL_H
#include <math.h>
#include <stdlib.h>
#include "mex.h"

/* OpenMP includes and compilation warnings */
#ifdef _OPENMP
#include <omp.h> /* Needed for VC++ linker */
#pragma message("Parallelism enabled in compilation!")
#else
#pragma message("Warning: Compile with the OpenMP flags to enable parallelism.")
#pragma message("\tSuggested command in GCC\n\t\tmex -v mutual.c CFLAGS=\"-fopenmp -Wall -O3 \\$CFLAGS\" LDFLAGS=\"-fopenmp \\$LDFLAGS\"")
#pragma message("\tSuggested command in VC++\n\t\tmex -v mutual.c COMPFLAGS=\"/openmp /O2 /Wall $COMPFLAGS\"")
#pragma message("\n\n")
#endif

#ifndef PI
#define PI 3.14159265358979323846
#endif
#define MU0 4*PI*1e-7
#define MUOVER4PI 1.0e-7
#define XX 0
#define YY 1
#define ZZ 2
#define EPS 1e-13
#define MAXsubfils 6          /* maximum subdivisions for subdividing filaments */

/* for deg_mutual, how small a dimension must be to be degenerate */
#define DEG_TOL 1.1e-4

/*-------------------------------------------------------------------------
 * Basic macro functions
 *------------------------------------------------------------------------*/

#define SQUARE(A) ((A)*(A))
#define CUBE(A) ((A)*(A)*(A))

#define nearzero(x) (fabs(x) < EPS)

#define TRUE (1==1)
#define FALSE (1==0)

#ifndef MAX
#define MAX(A,B)  ( (A) > (B) ? (A) : (B) )
#endif

#ifndef MIN
#define MIN(A,B)  ( (A) < (B) ? (A) : (B) )
#endif

#define compare(x,y,eps) (  (((x)==0 && (y)==0) || (fabs((x) - (y)) < eps*((x) + (y)) )) \
  ? 0 : (  ((x) > (y)) ? 1 : -1 )  )

#ifndef SIGN
#define SIGN(A)  ( (A) < 0.0 ? -1.0 : 1.0 )
#endif
#ifndef NORM
#define NORM(A)  ( (A) < 0.0 ? -(A) : (A) )
#endif
#ifndef PI
#define PI 3.14159265358979323846
#define PI2 1.570796326794897
#endif
  
/* These are missing in some math.h files */
#define atanh(x) 0.5*log((1+(x))/(1-(x)))
#define asinh(x) log((x)+sqrt((x)*(x)+1)) 
#define finite(x) ((x) != HUGE_VAL && (x) != -HUGE_VAL)
  
/*---------------------------------------------------------------------------------------- */
/* Utility functions (some of which were turned into macros) */
/*---------------------------------------------------------------------------------------- */
#define magdiff2(fil1, node1, fil2, node2) \
          ( SQUARE(fil1->x[node1] - fil2->x[node2]) \
	      +SQUARE(fil1->y[node1] - fil2->y[node2]) \
	      +SQUARE(fil1->z[node1] - fil2->z[node2]))
           
#define dotprod(fil1,fil2) \
          (  (fil1->x[1] - fil1->x[0])*(fil2->x[1] - fil2->x[0]) \
          + (fil1->y[1] - fil1->y[0])*(fil2->y[1] - fil2->y[0]) \
          + (fil1->z[1] - fil1->z[0])*(fil2->z[1] - fil2->z[0]) )
   
#define vdotp(v1,v2) \
          (v1[XX]*v2[XX] + v1[YY]*v2[YY] + v1[ZZ]*v2[ZZ])

#define dotp(x1,y1,z1,x2,y2,z2) \
          ((x1)*(x2) + (y1)*(y2) + (z1)*(z2))
     
#define getD(fil,D) \
          D[XX] = fil->x[1] - fil->x[0]; \
          D[YY] = fil->y[1] - fil->y[0]; \
          D[ZZ] = fil->z[1] - fil->z[0]
          
#define getr(x,y,z,s,t,D) \
          *x = s[XX] + t*D[XX]; \
          *y = s[YY] + t*D[YY]; \
          *z = s[ZZ] + t*D[ZZ]

#define mag(x1,y1,z1) \
          sqrt( (x1)*(x1) + (y1)*(y1) + (z1)*(z1) )

#define magsq(x1,y1,z1) \
          ( (x1)*(x1) + (y1)*(y1) + (z1)*(z1) )
          
#define fill_4(vec,E,a,d) \
          vec[0] = (E) - (a); \
          vec[1] = (E) + (d) - (a); \
          vec[2] = (E) + (d); \
          vec[3] = (E)

#define dist_between(x1,y1,z1,x2,y2,z2) \
          sqrt(SQUARE((x1) - (x2)) + SQUARE((y1) - (y2)) + SQUARE((z1) - (z2)))

#define aspectratio(fil) \
          ((fil->height >= fil->width) ? \
                  (fil->height/fil->width) :\
                  (fil->width/fil->height))  

typedef struct Filament {
  double x[2], y[2], z[2];  /* endpoints */
  double length, area, width, height;
  double *lenvect;        /* vector along the length of filament */
  double *widvect;
  double widthdir[3];
} FILAMENT; 

/* Counters
extern int num_mutualfil;
extern int num_perp;
extern int num_fourfil;
extern int num_exact_mutual;
extern int num_self;
extern int num_quadFil;*/

/* Filament utility functions */
extern void get_wid();
extern void get_height();
extern double min_endpt_sep();
extern double dist_betw_pt_and_fil();

/* Main working functions */
extern double mutual();
extern double self();
extern double cuber();
#endif /*MUTUAL_H*/