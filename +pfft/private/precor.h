#include <vector>
#include <algorithm> // for std::sort
#include <cstdlib> // for std::div
#include "mex.h"
#include "matrix.h" // for myAssert
#include "omp.h"

//-------------------------------------------------------------------------
// The following code tries to figure out whether a hash table is available 
// for the lookup table. 

// Visual Studio > 2010 or otherwise C++2011
#if (_MSC_VER > 1600) || (__cplusplus == 201103L)
#include <unordered_map>
#define MAP(TYPE) std::unordered_map<mwIndex,TYPE>
// Visual Studio 2003-2010
#elseif (_MSC_VER >= 1310 && _MSC_VER <= 1600) 
#include <hash_map>
#define MAP(TYPE) stdext::hash_map<mwIndex,TYPE>
#else
// Out of ideas
#include <map>
#define MAP(TYPE) std::map<mwIndex,TYPE> // Binary tree map!!
#endif
//-------------------------------------------------------------------------

// Gets the middle index
mwIndex getmid (const std::vector<mwIndex> &N) {
    mwIndex mid = 0;
    int cumprod = 1;
    for(int d=0; d<N.size(); d++) {
        cumprod *= N[d];
        mid += cumprod;
    }
    return (mid/2);
}

// Makes a list of subscripts according to raster indexing
std::vector<mwIndex> geti2i(const std::vector<mwIndex> &N) {
    // Count the total number of grid points + do some prep work
    std::vector<mwIndex> N2(N.size());
    std::vector<mwIndex> cumprod(N.size()+1,1); // for big grid
    mwIndex numel = 1; ldiv_t tmp;
    for (int d=0; d<N.size(); d++) {
        tmp = std::ldiv(N[d],2);
        
        // N is even. Then the grid is nicely segmented [0:N2-1,-N2:-1]
        if (tmp.rem == 0 )
            N2[d] = tmp.quot;
        else // N is odd. Throw an exception because our code can't handle it
            mexErrMsgTxt("The size of the Green's function matrix must be exactly twice that of the grid");
        
        numel *= N2[d];
        cumprod[d+1] = cumprod[d] * N[d]; // This is used for sub2ind conversion
    }
    
    // Iterate through each grid point.
    // Convert the index into subscripts then back to big grid indices
    std::vector<mwIndex> myi2i(numel,0); // init to zeros
    #pragma omp parallel default(shared)
    {
    ldiv_t tmp;
    std::vector<mwIndex> sub(N.size()+1); 
    #pragma omp for
    for(long int k=0; k<numel; k++) {
        sub[0] = k;
        for (int d=0; d<N.size(); d++) {
            // We use long division along each subscript, so that the method 
            // generalizes to multiple dimensions and can be highly 
            // parallelized
            tmp = std::ldiv(sub[d],N2[d]);
            
            // Obtain subscripts
            sub[d+1] = tmp.quot;
            sub[d] = tmp.rem;   
            
            // Add this component to index
            myi2i[k] += sub[d]*cumprod[d];
        }
    }
    } // #pragma omp parallel
    return myi2i;
}

// Stupid MATLAB mxAssert gets taken out when compiled without -g
void myAssert (bool val, const char *msg) {
    if (!val)
        mexErrMsgTxt(msg);
}