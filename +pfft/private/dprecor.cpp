/*
% This file is part of the matlab-pfft package by Richard Y. Zhang.
% http://web.mit.edu/ryz/www
% 
% Copyright (c) 2013, Richard Y. Zhang
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
% 1. Redistributions of source code must retain the above copyright notice, 
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright 
% notice, this list of conditions and the following disclaimer in the 
% documentation and/or other materials provided with the distribution.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR 
% A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
% HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
% TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * DPRECOR.CPP - Generalized Precorrection for real double matrices
 * function [s] = dprecor(G,P,I,i,j)
 *
 * This file aspires to be the most general precorrection routine possible, 
 * capable of processing an arbitrary number of dimensions, non-canonical 
 * projection / interpolation. 
 *
 *
 * The inputs are G,P,I,i,j and the matrix outputs a matrix s of the same 
 * size as i and j.
 * G -   Green's function block, of size 2M x 2N x 2O x 2P .....
 *       effectively, this is the block you do FFT onto.
 * P -   Projection matrix in sparse matrix form. The columns are the basis 
 *       functions and the rows are the grid points.
 * I -   TRANSPOSE of the interpolation matrix in sparse form. Again, the 
 *       columns are the testing functions and the rows are the grid points.
 * i,j - Column, row combinations to be precorrected. 
 *
 * This file is truly designed to be used with MATLAB. It works only with
 * MATLAB's sparse matrix format.
 *
 * Some efficiency is certainly tossed away for the sake of generality. By 
 * not forming a template 
 *
 * Below is the MATLAB-style pseudocode equivalent of the code below
 *
    function [s] = precor(G,P,I,i,j)  
    % Makes a list of subscripts sorted by their indices
    i2s = index_to_subscript_map(N); 

    for each i,j pair
        % Fetch their respective columns and row numbers
        [src_idx,~,thisP] = find(P(:,j)); % source
        [fld_idx,~,thisI] = find(I(:,i)); % field

        % Form convolution matrix for this particular pairing of grid points
        G = convmat(G,i2s(src_idx),i2s(fld_idx));
        
        % Perform projection and convolution to obtain grid potentials
        fld_pot = G*thisP;
        
        % Interpolate to obtain the precorrection
        s(ii) = thisI.'*fld_pot;
    end

    end
 */

#include "precor.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, 
  const mxArray *prhs[]) {
    myAssert(nrhs == 5 && nlhs == 1, "Incorrect number of inputs or outputs. function [s] = dprecor(G,P,I,i,j)");
    //---------------------------------------------------------------------
    // function [s] = precor(fftshift(G),P,I,i,j)
    // G - complex or real, double, full, Green's function
    // P - complex or real, double, sparse, projection matrix
    // I - complex or real, double, sparse, TRANSPOSE of interpolation matrix
    // i,j - real, double, sparse, interacting pairs
    // 
    // GREEN'S FUNCTION MATRIX
    myAssert(!mxIsComplex(prhs[0]), "The matrix G must be real. Call CPRECOR for processing complex matrices");
    myAssert(!mxIsSparse(prhs[0]), "The matrix G must be full");
    myAssert(mxIsDouble(prhs[0]), "The matrix G must contain doubles");
    double *G = mxGetPr(prhs[0]);
    const mwSize *tmp1 = mxGetDimensions(prhs[0]);
    mwSize tmp2 = mxGetNumberOfDimensions(prhs[0]);
    if ((tmp2 == 2) && (tmp1[1] == 1)) // Fix matlab's awkward dimensions counting
        tmp2 = 1;
    const std::vector<mwIndex> N(tmp1,tmp1+tmp2); 
    
    // Projection matrix
    // assert: size of P, P is real, P is double
    myAssert(!mxIsComplex(prhs[1]), "The matrix P must be real. Call CPRECOR for processing complex matrices");
    myAssert(mxIsSparse(prhs[1]), "The matrix P must be sparse");
    myAssert(mxIsDouble(prhs[1]), "The matrix P must contain doubles");
    double *prP = mxGetPr(prhs[1]);
    mwIndex *irP = mxGetIr(prhs[1]);
    mwIndex *jcP = mxGetJc(prhs[1]);

    // Interpolation matrix
    // assert: size of I, I is real, I is double
    myAssert(!mxIsComplex(prhs[2]), "The matrix I must be real. Call CPRECOR for processing complex matrices");
    myAssert(mxIsSparse(prhs[2]), "The matrix I must be sparse");
    myAssert(mxIsDouble(prhs[2]), "The matrix I must contain doubles");
    myAssert(mxGetM(prhs[2])==mxGetM(prhs[1]), 
            "The matrices P and I must contain the same number of rows. (I is the transpose of the interpolation matrix)");
    double *prI = mxGetPr(prhs[2]);
    mwIndex *irI = mxGetIr(prhs[2]);
    mwIndex *jcI = mxGetJc(prhs[2]);
    
    // Row and column indices
    // assert: number of elements is the same.
    mwIndex nument = mxGetNumberOfElements(prhs[3]);
    myAssert(nument == mxGetNumberOfElements(prhs[4]),
            "Number of elements in i and j must be the same");
    myAssert(!mxIsComplex(prhs[3]) && !mxIsComplex(prhs[4]),
            "The indices i and j must be real");
    double *i = mxGetPr(prhs[3]);
    double *j = mxGetPr(prhs[4]);
    
    // Prepare output
    plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[3]), mxGetN(prhs[3]),mxREAL);
    double *s = mxGetPr(plhs[0]);
    
    //---------------------------------------------------------------------
    // Generate a list to map index (small grid) -> index (big grid)
    // Use static language to prevent multiple initializations
    static std::vector<mwIndex> prevN;
    static std::vector<mwIndex> i2i; // List
    static mwIndex mid; // midpoint
    if (prevN != N) {
        i2i = geti2i(N);
        mid = getmid(N);
        prevN = N;
    }
    myAssert(i2i.size() == mxGetM(prhs[1]), 
            "The Green's function kernel provided must be the one used to perform the fft convolution");
    //---------------------------------------------------------------------
    
    #ifdef _OPENMP
        omp_set_num_threads(omp_get_num_procs());  
    #endif
    
    #pragma omp parallel default(shared)
    {
    // Private variable initialization
    mwIndex idx, src_bi, fld_ti, num_p, num_i, *src_idx, *fld_idx;
    double *thisP, *thisI, precor, fld_pot;
    
    // Optimization variables for each new column of P
    std::vector<std::pair<int,double> > pidxlst; // stores the i2i lookup, and projection value
    MAP(double) pot_map; // Stores potential evaluated on the grid
    
    #pragma omp for
    for (long int k=0; k<nument; k++) {
        src_bi = (mwIndex) j[k] - 1; // Basis fn id
        fld_ti = (mwIndex) i[k] - 1; // Testing fn id
        // To do: assert that these numbers are within range
        
        // The first entries for these columns
        mwIndex src_k = jcP[src_bi];
        mwIndex fld_k = jcI[fld_ti];
        
        // Shift pointers to the beginning of these columns
        thisP = &prP[src_k]; // Entry
        thisI = &prI[fld_k];
        src_idx = &irP[src_k]; // Row index (i.e. Grid index)
        fld_idx = &irI[fld_k];
        
        // Figure out how many elements are in these columns
        num_p = jcP[src_bi+1] - src_k;
        num_i = jcI[fld_ti+1] - fld_k;
        
        // Prepare the optimization lookup table and vector
        if (pidxlst.empty() || j[k] != j[k-1]) { // "if processing new column of P"
            // Clear lookup table
            pot_map.clear();
            
            // Rewrite over the vector, do not bother to clear
            pidxlst.reserve(num_p); // Allocate enough resources to not bug out.
            for(mwIndex j=0; j<num_p; j++) {
				// pidxlst[j].first = index of the "from"
				// pidxlst[j].second = value of the projection point
                pidxlst[j].first = i2i[src_idx[j]] - mid;
                pidxlst[j].second = thisP[j];
            }
            // Sort to make it slightly faster to iterate through
            std::sort (pidxlst.begin(), pidxlst.begin()+num_p);
        }
        
        // Form Precorrection    
        // Iterate through each affected grid point
        precor =0;
        for(mwIndex i=0; i<num_i; i++) { // each field point
            // fld_pot is the potential evaluated at this grid point
            fld_pot = pot_map[fld_idx[i]]; // Look for it in our table
            if (fld_pot == NULL) { // Not found
                mwIndex this_fld_idx = fld_idx[i];
                // Find where it is on the big grid
				mwIndex new_fld_idx = i2i[this_fld_idx]; 
                // Compute by convolution dot product
                // (Yes, randomly accessing G. It's somewhat ordered due 
                // to the presorting)
                fld_pot = 0;
                for(mwIndex j=0; j<num_p; j++)
                    fld_pot += G[new_fld_idx - pidxlst[j].first]
                                  * pidxlst[j].second;
                pot_map[this_fld_idx] = fld_pot; // Insert into the table
            }
            // Interpolation inner product
            precor += thisI[i] * fld_pot;
        }
        s[k] = precor; // store
    }
    } // #pragma omp parallel
}