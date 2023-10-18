% This file is part of the matlab-pfft package by Richard Y. Zhang.
% http://web.mit.edu/ryz/www
% 
%PFFT MATLAB Precorrected FFT version 2013
% Author: Richard Y. Zhang (ryz@mit.edu)

% Consider the n-body problem, for a set of points in 3-space:
%             x y z
%   points = [1 2 3;
%              ...
%             3 3 2];
% each weighed according to some mass:
%        x = [1, 2, ... 5]';
% The standard method to computing the potential at all of these points is 
% to form the Green's coupling matrix A, where, looping over the rows i and 
% columns j:
%   A(i,j) = 1/norm(points(i,:) - points(j,:));
% then performing the matrix-vector product:
%   >> b = A * x; 
% to yield b as the vector of potentials at each of the points.

% In fast solvers, this same calculation is accelerated. First we perform a
% setup, with 1/r as the Green's function:
%   >> p = pfft(@inv_r, points);
% and then perform the mv-product approximately:
%   >> b = p.fastmv(x);

classdef pfft < handle
    properties
        % Grid
        pfft_grid;
                
        % Kernel dependent terms
        Tker; % Toeplitz kernel tensor
        Hker; % Hankel kernel tensor
        rotidx = {}; % Hankel rotation indices
        outidx = {}; % Output isolation indices
        
        % Matrices
        Pmat; % Projection Matrix
        Imat; % Interpolation Matrix
        PCmat; % Precorrection Matrix
        sym = false; % Symmetry allows some savings in computation and memory
        
        % PFFT OPTIONS
        ord=2; % Polynomial expansion order. Higher is more accurate.
        points_per_grid=1; % Number of grid points to allocate per quadrature point
        basis_quad_ratio=1; % Number of quadrature points required to increase basis count by 1
        batchsize = 10e6; % Size of each batch of precorrection / direct stencil computation
        gpu=false; % Make use of gpuArray capability? (Only available on CUDA machines)
        verbose=true; % mode 0 for clean, mode 1 for verbose. 
        
        % Internal warnings
        dsetup=false; % Has the direct matrix been set up?
    end
    
methods
%--------------------------------------------------------------------------
% PFFT Constructor - Sets up a pFFT object
%--------------------------------------------------------------------------
    function hObj = pfft(green, testing, basis, setdirect)
        if nargin == 2 || isempty(basis)
            hObj.sym = true;
            basis = [];
        end
        if nargin < 4
            setdirect = true;
        end
        hObj.msgtime; % Start timer
        tparticle = isnumeric(testing); bparticle = isnumeric(basis);
        
        %------------------------------------------------------------------
        % GRID SETUP
        %------------------------------------------------------------------       
        % TODO: Wavelength restricted grid, for e.g. eikr kernel
        hObj.optim_grid(testing,basis,hObj.ord); % Generate the optimal grid.
        hObj.set_kernel(green); % Set up the Toeplitz kernel, to be precorrected.
        
        %------------------------------------------------------------------
        % PERFORM THE PROJECTION
        %------------------------------------------------------------------
        % Project and interpolate and precorrect
        I = hObj.calc_proj(testing);
        if ~hObj.sym
            P = hObj.calc_proj(basis);
            if setdirect
                hObj.set_geom(P,I');
            else
                hObj.Pmat = P;
                hObj.Imat = I';
            end
        else
            if setdirect
                hObj.set_geom(I);
            else
                hObj.Pmat = I;
            end
        end
        
        %------------------------------------------------------------------
        % Compute nearby interactions using the Green's function provided.
        % (This is only applicable to particle-particle problems.)
        %------------------------------------------------------------------
        if bparticle && tparticle
            if ~hObj.sym
                dirfun = @(i,j) green(num2cell(testing(i,:)-basis(j,:),1));
            else
                dirfun = @(i,j) green(num2cell(testing(i,:)-testing(j,:),1));
            end
            hObj.calc_direct(dirfun);
        end
    end
%--------------------------------------------------------------------------
% Grid Optimizer - sets up the optimal grid in order to minimize costs.
%--------------------------------------------------------------------------
    function optim_grid(hO, testing, basis, ord)
        hO.msgtime;
        
        % Extract bounds information for the testing function
        
        if isnumeric(testing)
            ubnd = max(testing,[],1);
            lbnd = min(testing,[],1);
            numt = size(testing,1);
            dims = size(testing,2);
        elseif isstruct(testing) && isfield(testing,'points') && isfield(testing,'weight')
            dims = size(testing(1).points,2);
            ubnd = -inf(1,dims); lbnd = inf(1,dims);
            for wi = 1:numel(testing)
                ubnd = max([testing(wi).points;ubnd],[],1);
                lbnd = min([testing(wi).points;lbnd],[],1);
            end
            numt = size(testing(1).points,1);
        else
            error('Wrong input format for testing function!');
        end
        
        % Continue expanding the box for the basis function
        if isnumeric(basis)
            ubnd = max([basis;ubnd],[],1);
            lbnd = min([basis;lbnd],[],1);
            numb = size(basis,1);
        elseif isstruct(basis) && isfield(basis,'points') && isfield(basis,'weight')
            for wi = 1:numel(basis)
                ubnd = max([basis(wi).points;ubnd],[],1);
                lbnd = min([basis(wi).points;lbnd],[],1);
            end
            numb = size(basis(1).points,1);
        elseif ~isempty(basis)
            error('Wrong input format for basis function!');
        end
        
        % Sample 1e5 random basis pairs and compute their L_inf distance.
        % Obtain the cumulative distribtion function
        [dist, cdf] = dist_cdf(testing, basis, 1e5);
        
        % grid_dist * padding * 2 is the max distance between neighbors
        gdist = dist(:) / ord;
        % Now cdf(gdist) gives number of number of neighbors for a
        % particular grid size
        
        % To compute fft costs, see pgrid. The number of grid points for 
        % a given grid delta is,
        ext = ubnd - lbnd;
        N = ceil(bsxfun(@rdivide,ext,gdist(:)))+1+ord/2;
        % and the big grid is twice this size in each dimension,
        N = prod(N*2,2);
        % The cost of FFT is ~ 4NlogN - 6N + 8 for split radix
        cfft = 4 .* N .* log2(N) - 6.*N + 8;
        
        % Optimize only for fast mv products.
        ctot = 100*cdf + cfft;
        [cmin,imin] = min(ctot);
        d0 = gdist(imin);
        
%         if hO.verbose
%             loglog(gdist,cdf,gdist,cfft,gdist,ctot);
%             legend('Direct stencil cost','FFT cost','Combined cost');
%             ylabel('# multiplications and additions');
%             xlabel('Grid spacing');
%             title(sprintf('Minimum cost: %d with grid delta: %g',cmin,d0));
%             pause
%         end
        
        % Output the grid
        hO.pfft_grid = pfft.pgrid(d0,[ubnd;lbnd],ord/2);
        numnei = hO.pfft_grid.nei(testing,basis);
        
        % Report the stats
        hO.msgtime('Grid prep'); % Display message
        hO.msg(['Grid size: ' mat2str(hO.pfft_grid.N) ' with '...
            mat2str(hO.pfft_grid.numel) ' points']);
        hO.msg(sprintf('Number of neighbors: %d',numnei));
        
    end
%--------------------------------------------------------------------------
% Convolution Kernel Setup - sets up the Toeplitz and Hankel kernels upon 
% the grid defined above.
%
% Every time this function is called, it re-sets-up all the Green's
% functions. This can be useful for cases where:
% 1. Only the Green's function needs to be updated, for example, in 
%    frequency-dependent sweeps. 
% 2. Precorrection needs to be performed using a different Green's function
%    to the ones used in the convolutions
%
% t_fun is the Toeplitz (r-r') function. h_fun is the Hankel (r+r')
% function. rotdims can be used to define a kernel that is partially
% hankel, e.g. (x-x',y-y',z+z'), let rotdims = [0 0 1];
%--------------------------------------------------------------------------
    function set_kernel(hO, t_fun, h_fun, rotdims)
        hO.msgtime; % Start timer
        
        % Delete the previous kernels stored
        hO.Hker = []; hO.Tker = [];
        if nargin == 1 
            hO.msg('Convolution kernels have been cleared.');
            return;
        end
        
        % Import variables
        g = hO.pfft_grid;
        
        % Convolution setup. This sets up the relevant convolution kernels.
        if nargin > 1 && ~isempty(t_fun)
            toeplitz_ker = t_fun(g.getCoordBox); 
            
            % remove non-finite values
            nfin_id = ~isfinite(toeplitz_ker);
            if any(nfin_id(:))
                warning('Non-finite quantities encountered in toeplitz kernel. They have been set to zero, but results may be wrong.');
                toeplitz_ker(nfin_id) = 0;
            end
            
            hO.Tker = fftn(toeplitz_ker);
            if hO.gpu % gpu option enabled?
                hO.Tker = gpuArray(hO.Tker);
            end
        end
        if nargin > 2 && ~isempty(h_fun)
            hankel_ker = h_fun(g.getCoordBox(rotdims));

            % remove non-finite values
            nfin_id = ~isfinite(hankel_ker);
            if any(nfin_id(:))
                warning('Non-finite quantities encountered in hankel kernel. They have been set to zero, but results may be wrong.');
                hankel_ker(nfin_id) = 0;
            end
            
            % Create the hankel rotation matrix
            hO.rotidx = repmat({':'}, 1, g.dims);
            for this_dim = 1:g.dims
                if rotdims(this_dim)
                    hO.rotidx{this_dim} = g.N(this_dim):-1:1;
                end
            end

            hO.Hker = fftn(hankel_ker);
            if hO.gpu % gpu option enabled?
                hO.Hker = gpuArray(hO.Hker);
            end
        end
        
        % Create selection index block. This is used during fastmv, 
        % to isolate the half / quarter / eighth of the toeplitz product
        % from the circulant product.
        for this_dim = 1:g.dims
            hO.outidx{this_dim} = 1:g.N(this_dim);
        end

        % Display message
        hO.msgtime('Kernel prep');
    end
%--------------------------------------------------------------------------
% Projection & Interpolation Setup - copies the projection & interpolation
% matrices, and performs precorrection on them using the currently stored
% Toeplitz Green's function.

% Since the precorrection is always unconditionally tied to the P & I
% matrices, we will take this opportunity to precorrect
% TODO: memory reduction for symmetric P&I cases
%--------------------------------------------------------------------------
    function PCmat = set_geom(hO, P, I, ij)
        if nargin == 4 % Allow grabbing precor
            if isempty(P)
                P = hO.Pmat;
            end
            if isempty(I) && ~hO.sym
                I = hO.Imat;
            end
        end
        
        assert(issparse(P),'P must a sparse matrix!');
        hO.Pmat = P;
        if nargin > 2 && ~isempty(I)
            assert(issparse(I),'I must a sparse matrix!');
            assert(size(P,1) == size(I,2), 'Number of grid points must be the same, i.e., P must have the same number of rows as I does columns.');
            hO.Imat = I;
        end
        
        % Skip setting up the precorrection if Toeplitz kernel is empty
        if isempty(hO.Tker) 
           hO.msg('No precorrection performed because Toeplitz kernel is empty.');
           return; 
        end
        
        % Start timer
        hO.msgtime;

        % Get the kernel back, rotate so that the middle is in the center
        G = fftshift(ifftn(hO.Tker));
        if hO.gpu
            G = gather(G); end
        
        if hO.sym
            % Figure out where to precorrect;
            if nargin <4
                [di,dj] = find(P.'*P);           
            else
                di = ij(:,1); dj = ij(:,2); assert(size(ij,2)==2);
            end
            hO.msg(['Number of neighbors: ' mat2str(numel(di))]);

            % Call the external mex functions to perform the precorrection
            if 1e8*norm(imag(G(:))) < norm(real(G(:)))
                pcfun = @(i,j) dprecor(real(G),P,P,i,j); % real double routine
            else
                pcfun = @(i,j) cprecor(G,P,P,i,j); % complex double routine
            end
            
            % Perform batch computation
            [pc,di,dj] = hO.batchcalc(pcfun,di,dj);
            
            % Chose output
            if nargout == 0
                hO.PCmat = sparse(di,dj,pc,size(P,2),size(P,2));
            else
                PCmat = sparse(di,dj,pc,size(P,2),size(P,2));
            end
        else % Non-symmetric
            % Figure out where to precorrect;
            if nargin <4
                [di,dj] = find(I*P);           
            else
                di = ij(:,1); dj = ij(:,2); assert(size(ij,2)==2);
            end
            hO.msg(['Number of neighbors: ' mat2str(numel(di))]);

            % Call the external mex functions to perform the precorrection
            if 1e8*norm(imag(G(:))) < norm(real(G(:)))
                pcfun = @(i,j) dprecor(real(G),P,I.',i,j); % real double routine
            else
                pcfun = @(i,j) cprecor(G,P,I.',i,j); % complex double routine
            end
            
            % Perform batch computation
            [pc,di,dj] = hO.batchcalc(pcfun,di,dj);
            
            % Chose output
            if nargout == 0
                hO.PCmat = sparse(di,dj,pc,size(I,1),size(P,2));
            else
                PCmat = sparse(di,dj,pc,size(I,1),size(P,2));
            end
        end
        
        % Display message
        hO.msgtime('Precorrection');
    end
%--------------------------------------------------------------------------
% Computes the projection of basis or points onto the pfft grid 
%--------------------------------------------------------------------------
    function Pmat = calc_proj(hO, basis)
        hO.msgtime;
        Pmat = hO.pfft_grid.proj(basis);
        hO.msgtime('Projection'); % Display timing message
    end
%--------------------------------------------------------------------------
% Direct Matrix Setup
% input: A function that computes elements of the A matrix
%--------------------------------------------------------------------------
    function Dmat = calc_direct(hO, Afun, ij)
        % TO DO: Break the computation into blocks for many computations
        if nargin <3 && isempty(hO.PCmat)
            error('Cannot compute direct interactions before precorrection');
        end
        hO.msgtime;
        if nargin == 3
            di = ij(:,1); dj = ij(:,2); assert(size(ij,2)==2);
        else
            [di,dj] = find(hO.PCmat);
        end
        
        % Perform batch calculation
        [ad,di,dj] = hO.batchcalc(Afun,di,dj);
        if hO.sym
            Dmat = sparse(di,dj,ad,size(hO.Pmat,2),size(hO.Pmat,2));
        else
            Dmat = sparse(di,dj,ad,size(hO.Imat,1),size(hO.Pmat,2));
        end
        
        % Combine the direct matrix with the precorrection if the user does
        % not collect it.
        if nargout == 0
            hO.msg('Direct Matrix combined with Precorrection');
            hO.PCmat = hO.PCmat - Dmat;
        else
            hO.msg('Direct Matrix collected by user');
        end
        hO.dsetup = true; % The user has set up a valid direct matrix.
        hO.msgtime('Direct matrix');
    end
%--------------------------------------------------------------------------
% Fast Matrix-Vector Product
%--------------------------------------------------------------------------
    function [mvprod, potfld] = fastmv(hO,x)
        % Error checks & data extraction
        if ~hO.dsetup % direct matrix
            warning('The nearby interactions were never computed. These results are probably wrong.'); end
        if isempty(hO.Pmat) || (isempty(hO.Imat) && ~hO.sym) % projection / interpolation
            error('The Projection / Interpolation matrices were never set up properly.'); end
        N = hO.pfft_grid.N; % Size of the kernel
        dims = hO.pfft_grid.dims;
        
        % Project
        Q = hO.Pmat*x; % Project onto pfft grid
        if hO.gpu % Move the projection to the GPU
            Q = gpuArray(Q); end
        if dims > 1
            Q = reshape(Q,N); end
        
        % Convolve in the Fourier Domain
        potfld = 0;
        if ~isempty(hO.Tker)
            potfld = potfld + ...
                fftn(Q,size(hO.Tker)) .* hO.Tker; 
        end
        if ~isempty(hO.Hker)
            potfld = potfld + ...
                fftn(Q(hO.rotidx{:}),size(hO.Hker)) .* hO.Hker; 
        end
        
        % Invert back to spatial domain
        potfld = ifftn(potfld);
        potfld = potfld(hO.outidx{:}); % Isolate the one corner
        if hO.gpu % Collect the product back from the gpu
            potfld = gather(potfld); end
        
        % Interpolate, resolve direct terms
        if hO.sym
            mvprod = hO.Pmat.'*potfld(:);
        else
            mvprod = hO.Imat*potfld(:);
        end
        if ~isempty(hO.PCmat) % Resolve the direct matrix
            mvprod = mvprod - hO.PCmat*x; end
    end
%--------------------------------------------------------------------------
% Timer + Message Display (Only works when in verbose mode)
%--------------------------------------------------------------------------
    function msgtime(hO,prompt)
        persistent t_verb
        if nargin == 2
            hO.msg(sprintf('%s: \t%g',prompt,toc(t_verb)));
        end
        t_verb = tic;
    end
    function msg(hO,prompt)
        if hO.verbose == 1
            fprintf('%s\n',prompt);
        end
    end
%--------------------------------------------------------------------------
% Execute fun(i,j) where i and j are long column vectors in batches.
% If hO.sym flag is true, then we know that: 
%  -- for every (i,j) pair there is a (j,i) pair somewhere
%  -- fun(i,j) = fun(j,i)
% This allows us to save some @fun executions.
%--------------------------------------------------------------------------
function [out,di,dj] = batchcalc(hO,fun,i,j)
        if hO.batchsize < 1e3
            error('Batch size defined too low!'); end
        assert(numel(i) == numel(j), 'Number of froms and tos must be the same');
        
        % Isolate the upper triangular pairs if symmetric.
        if hO.sym
            % Find pairs
            uppertriag = (i > j);
            selfterms = (i == j);
            
            % Make new i,j pairs
            slf = i(selfterms);
            di = i(uppertriag); dj = j(uppertriag);
        else
            di = i; dj = j;
        end
            
        %------------------------------------------------------------------
        % Do the first x number in batches
        numbatch = floor(numel(di)/hO.batchsize);
        numtodo = numbatch * hO.batchsize;
        
        % Put the first numtodo in a rect matrix
        % Chop these off from i and j
        bi = reshape(di(1:numtodo),hO.batchsize,numbatch); %"batch"
        bj = reshape(dj(1:numtodo),hO.batchsize,numbatch);
        ri = di(numtodo+1:end); rj = dj(numtodo+1:end); %"remainder"
        
        % Begin by filling out the ones allocated to batches
        out = zeros(size(bi)); fprintf('\n');
        for ii = 1:numbatch % You can change this to parfor
            out(:,ii) = fun(bi(:,ii),bj(:,ii));
            fprintf('*');
            if mod(ii,20) == 0
                fprintf('\n'); end
        end
        fprintf('\n');
        out = out(:);
        
        % Tag on the last batch
        if ~isempty(ri)
            out = [out;fun(ri,rj)];
        end
        %------------------------------------------------------------------
        
        % Copy upper triag to lower triag if symmetric
        if hO.sym
            % Pass message
            hO.msg('Utilized symmetry to reduce computations');
            
            selfout = fun(slf,slf);
            di = [di(:),dj(:)]; dj = di(:,[2,1]);
            di = [di(:);slf]; dj = [dj(:);slf];
            out = [out;out;selfout];
        end
    end
end
end
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