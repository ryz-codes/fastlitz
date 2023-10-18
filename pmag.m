% This file is part of the pmag package by Richard Y. Zhang.
% http://web.mit.edu/ryz/www
% 
% This is the object version. It is designed so that setup can be
% performed on a fast multicore computer, and individual extraction can be
% done on a slow computer.
classdef pmag <handle
    properties
        % Geometry properties
        geom; % Problem geometry
        Lvec; %Length vectors. Used in MV products.
        Nfils; %Number of filaments per segment
        Nsegs; %Number of segments per terminal
        Nterms; %Number of terminals
        
        % Core Matrix-Vector Product variables
        p; %pfft object
        Dmat; %Direct matrix
        Rmat; %Resistance matrix
        Lblk; Rblk; %Preconditioners
        
        % Extraction variables
        M; %Mesh Matrix. The first Nterms rows correspond to terminals.
        S; % Terminal matrix
        Marea; %Mesh area vector, used for B.da computations
        Mcen; %Centroids of the branches
        
        % Options
        ord=2; %PFFT projection / interpolation order
        verbose=true; %Printouts
    end
methods
function [K] = pmag(geom)
    if nargin ==0
        return; % Default constructor
    end
    O = geom.O; L = geom.L; W = geom.W; H = geom.H;

    % Check the problem type - multiterminal or single terminal
    % TODO: Add more comprehensive error checks
    if (iscell(O) && iscell(L) && iscell(W) && iscell(H))
        if (iscell(O{1}) && iscell(L{1}) && iscell(W{1}) && iscell(H{1}))
            % do nothing
        else
            error('Single-terminal extraction deprecated');
        end
    else
        error('No-terminal extraction not implemented');
    end

    % Expand out to branches, generate mesh matrices
    K.Nterms = numel(O);
    K.Nfils = mean(cellfun(@(x)size(x,1),O{1}));
    K.Nsegs = zeros(K.Nterms,1);
    for ij = 1:K.Nterms
        tmpchk = cellfun(@(x)size(x,1),O{ij});
        assert(all(K.Nfils == tmpchk),'Must have same number of filament in each segment!!')
        K.Nsegs(ij) = numel(O{ij});

        % Expand out the branches
        O{ij} = vertcat(O{ij}{:}); L{ij} = vertcat(L{ij}{:}); 
        W{ij} = vertcat(W{ij}{:}); H{ij} = vertcat(H{ij}{:});
    end
    O = vertcat(O{:}); L = vertcat(L{:}); 
    W = vertcat(W{:}); H = vertcat(H{:});
    
    % Recover path variables, for pFFT and mesh areas
    K.Lvec = L;
    K.Mcen = O+W/2+H/2;
    
    % Mesh matrix, assign terminals
    MT = termMat(K.Nsegs,K.Nfils);
    K.M = [MT;loopMat(K.Nsegs,K.Nfils)];
    K.S = false(K.Nterms + (K.Nfils-1)*sum(K.Nsegs),1);
    K.S(1:K.Nterms) = true; % Terminal matrix
    
    % Total numbers
    Nbranch = size(O,1);
    
    %----------------------------------------------------------------------
    % Compute proximity effect area vector.
    %----------------------------------------------------------------------
    seg = 1:K.Nfils:Nbranch;
    
    % Strand level loops
    bfwd = repmat(seg,K.Nfils-1,1); % Loop heads
    brev = bfwd + repmat((1:K.Nfils-1).',1,sum(K.Nsegs)); % return paths
    sega = surfint(bfwd,brev); % Surface integral for strand-level loops
    
    % Iterate through terminal heads
    tla = zeros(K.Nterms,3); % Terminal LOOP voltages
    for ti = 2:K.Nterms
        % Loop the first term with the rest
        sfwd = find(MT(1,:));
        srev = find(MT(ti,:));
        % Sum over each pair of filaments.
        tla(ti,:) = sum(surfint(sfwd,srev));
    end
    T = [termMat(K.Nterms,1);loopMat(1,K.Nterms)];
    terma = T\tla; % Individual terminal voltages
    K.Marea = [terma;sega];

    % Diagonal Resistance matrix
    K.Rmat = K.resist(ones(Nbranch,1),[],L,W,H);

    % Set up the preconditioner
    [Lself] = K.induct(O,L,W,H); % Precorrect, get preconditioner
    K.Lblk = Lself(1:K.Nfils,1:K.Nfils); % Same self block!
    K.Rblk = diag(K.Rmat(1:K.Nfils));
    
    % Display
    fprintf('Setup complete for %d elements...\n',Nbranch);
    
    
    function a = surfint(fwd,rev)
        fwd = fwd(:); rev = rev(:);
        % Computes the area vector for each pair of fwd & rev branches.
        % The fwd and rev are given by their branch indices.
        w_ = K.Mcen(rev,:) - K.Mcen(fwd,:);
        p_ = w_ + K.Lvec(rev,:);
        q_ = w_ - K.Lvec(fwd,:);
        a = 0.5 * cross(p_,q_,2);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Impedance extraction
% Options: {'Ztot'} | 'Zmat' | val
%   'Ztot' - Total impedance 
%   'Zmat' - Impedance matrix for all terminals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Zt,s,Ib] = impedance(K,f,opts)
    w = 2*pi*f;
    
    % Solve for each frequency
    Zt = zeros(numel(K),numel(w));
    s = cell(numel(K),numel(w));
    Ib = cell(numel(K),numel(w));
    
    % Sweep the models
    for i = 1:numel(K)
    % Objective vector
    Vm0 = zeros(size(K(i).M,1),1); 
    Vm0(K(i).S) = 1; % Set all terminals to 1V.
    
    % Sweep the frequencies    
    Zti = zeros(1,numel(w));
    si = cell(1,numel(w));
    Ibi = cell(1,numel(w));
    for j = 1:numel(w)
        if K(i).verbose
            fprintf('\nStartin f = %e\n',f(j));
        end
        
        % Extract the impedance
        [Im0,si{j},Ib0] = K(i).solve(w(j),Vm0);
        Zti(j)=1/sum(Im0(K(i).S)); 
        
        % Store currents
        Ibi{j} = Ib0;
        
        if K(i).verbose
        fprintf('\nIterations: %g, Residual: %e, Time: %d\n',...
            si{j}.iters,si{j}.resid,si{j}.time);
        end
    end
    
    % Save
    Zt(i,:) = Zti; s(i,:) = si; Ib(i,:) = Ibi;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Power extraction
% Options: {'T'} | 'A/m'
% Input can be a row vector or a function that takes row vectors of (x,y,z)
% and output row vectors of (Bx,By,Bz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P s Ib] = power(K,f,exc,opts)
    w = 2*pi*f;
    
    % Solve for each frequency
    P = zeros(numel(K),numel(w));
    s = cell(numel(K),numel(w));
    Ib = cell(numel(K),numel(w));
    
    for i = 1:numel(K)
    % Objective vector
    % Verify B field data
    if isnumeric(exc) && numel(exc) ==3
        texc = exc(:).';
    elseif isa(exc, 'function_handle')
        [Bx,By,Bz] = exc(K(i).Mcen(:,1),K(i).Mcen(:,2),K(i).Mcen(:,3));
        texc = [Bx,By,Bz]; % Evaluate field function at loop centroids
    end
        
    % dot product B with area
    Vm0 = -sum(bsxfun(@times,K(i).Marea,texc),2); 
    
    % Sweep the frequencies
    Pi = zeros(1,numel(w));
    si = cell(1,numel(w));
    Ibi = cell(1,numel(w));
    for j = 1:numel(w)
        if K(i).verbose
        fprintf('\nStartin f = %e\n',f(j));
        end
        % Extract the impedance
        [~,si{j},Ib0] = K(i).solve(w(j),1j*w(j)*Vm0);
        % Compute the power dissipated in each branch
        % Pb = 0.5 * |Ib|^2 * Rb (sinusoids)
        % Pout = sum(Pb)
        % Rac = Pout / Iout^2
        Pi(j) = sum(0.5*K(i).Rmat.*abs(Ib0).^2);
        
        % Store currents
        Ibi{j} = Ib0;
        
        if K(i).verbose
        fprintf('\nIterations: %g, Residual: %e, Time: %d\n',...
            si{j}.iters,si{j}.resid,si{j}.time);
        end
    end
    
    % Store
    P(i,:) = Pi; s(i,:)=si; Ib(i,:) = Ibi;
    
    end
    % Scale the excitation to H
    if nargin==4 && isequal(opts,'A/m')
        mu0 = 4*pi*1e-7;
        P = P *mu0^2;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTATIONAL STUFF BELOW, DO NOT TOUCH
% Preconditioner: Use left-preconditioner instead of right.
% Given inv(P)Ax = inv(P)b, we set c = inv(P)b by solving Pc = b.
% Then, at each iteration, we have inv(P)Ax1 = c1. We do Mv production
% first to give Ax1 = b1, then we solve to get c1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Im0, stats, Ib0] = solve(K,w,Vm0)  
    % Solve iteratively
    st = tic;
    PRECON = true;
    
    % Compute the preconditioner diagonal
    Zself = K.Rblk+1j*w*K.Lblk;
    Yself = inv(Zself);
    
    % Solve by GMRES
    if PRECON % Preconditioner is MZselfM'
        % Iteratively solve the decoupled 
        pc = @(x) K.M*blkmult(Zself,(K.M'*x));
        [Vm0,~] = gmres(pc,Vm0,[],1e-6,50);
        Vm0 = full(Vm0);
    end
    [Im0, ~,residual,it] = gmres(@PMZMT,Vm0,[],1e-6,200);
    
    % Solve for branch currents.
    Ib0 = K.M'*Im0;
    
    % Collect statistics
    stats.time = toc(st);
    stats.iters = max(it);
    stats.resid = residual(1);

    % Multiply by mesh impedance matrix, right preconditioner
    function Vm = MZMTP(Im)
        if PRECON
            % Do preconditioner
            Ib = blkmult(Yself,K.M\Im);
        else
            Ib = K.M'*Im;
        end
        if w > 0 % Don't bother with pfft if you're just gonna do DC.
            Vm = K.M*(1j*w*K.induct(Ib)+K.Rmat.*Ib);
        else
            Vm = K.M*(K.Rmat.*Ib);
        end
        fprintf('*');
    end
    
    % Multiply by mesh impedance matrix, left preconditioner
    function Vm = PMZMT(Im)
        Ib = full(K.M'*Im);
        if w > 0 % Don't bother with pfft if you're just gonna do DC.
            Vb = (1j*w*K.induct(Ib)+K.Rmat.*Ib);
        else
            Vb = (K.Rmat.*Ib);
        end
        if PRECON
            % Do preconditioner
            Ib = blkmult(Yself,Vb);
            % Hack below, back substituting M' was very slow. To do:
            % replace with precomputed LU factorization.
            [Vm,~] = lsqr(K.M',Ib,1e-6,200);
        else
            Vm = K.M*Vb;
        end
        fprintf('*');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [V] = induct(K,varargin)
%INDUCT: calculates the PEEC matrix vector product [L]*I for a flat
% coil systems, using the pre-corrected FFT algorithm.
%   
%   Usage: 
%   1.  For each new set of filaments, to initialize the solver grid kernel, 
%       and precorrection by calling
%           induct(O,L,W,H);
%       a grid is then constructed for approximately 1 grid cell per 3
%       elements. Precorrection is performed for these filament sizes.
%
%   1a. If the direct matrix is desired, for example as a preconditioner,
%       then call INSTEAD:
%           D = induct(O,L,W,H);
%
%   2.  To perform MV products, do the following:
%           V = induct(I).
%   
%   2a. If a set of subtracted kernels are to be added, then call
%           induct(g); 
%       where g is the subtracted kernel to be added. WARNING: g must not
%       include the 1/r portion (it will be added automatically). This will 
%       leave the precorrection unaffected.
% Internal settings

% Get the PFFT Package
import pfft.*

% Initialization routine for free-space elements
if (nargin==5)
    % Process Geometry
    O = varargin{1};
    L = varargin{2};
    W = varargin{3};
    H = varargin{4};
    
    % Setup basis function quadrature points using sparse grid quadrature
    [nodes, qwei] = nwspgr('KPU', 3, K.ord);
    Nq = size(nodes,1);
    basis(Nq) = struct;
    for ii = 1:Nq
        basis(ii).points = O + bsxfun(@times,L,nodes(ii,:)) + bsxfun(@times,W,nodes(ii,:)) ...
                    + bsxfun(@times,H,nodes(ii,:));
        basis(ii).weight = qwei(ii);
    end
    
    % Setup pfft
    K.p = pfft(@inv_r,basis);
    clear basis nodes qwei;
    
    % Compute the direct matrix, collect externally.
    Afun = @(to,from) 1e7*mutual([O(from,:),O(to,:)]',...
        [L(from,:),L(to,:)]',...
        [W(from,:),W(to,:)]',...
        [H(from,:),H(to,:)]');
    K.Dmat = K.p.calc_direct(Afun);
    
    % Export the preconditioner
    V = K.Dmat*1e-7;
    
    % Steal PC mat
    % Sanity check. nnz(Dmat) must be same as nnz(prec)
    prec = K.p.PCmat; K.p.PCmat = [];
    [di,dj,dk] = find(K.Dmat); [pi,pj,pk] = find(prec);
    assert(isequal(di,pi) && isequal(dj,pj));
    
    % Merge PC with Dmat
    for dim = 1:3
        dk = dk - L(pi,dim).*pk.*L(pj,dim);
    end
    K.Dmat = sparse(di,dj,dk,size(K.Dmat,1),size(K.Dmat,2));
    
% Update the subtracted kernel
elseif nargin == 2 && nargout == 0
    tph_obj = varargin{1};
    assert(isa(tph_obj,'tph'),'Must supply a tph object!');
    
    % Setup green's functions
    t_fun = @(in) 1e7*tph_wrap(tph_obj,'T',in) + inv_r(in);
    h_fun = @(in) 1e7*tph_wrap(tph_obj,'H',in);

    % Init the new TPH kernel
    K.p.init_kernel(t_fun, h_fun, [0 0 1]);
    
% Single MV product
elseif nargin == 2 && nargout == 1
    V = 0;
    for ii = 1:3 % Perform pFFT dot product
        V = V + K.Lvec(:,ii) .* (K.p.fastmv(varargin{1}.*K.Lvec(:,ii)));
    end
    V = V + K.Dmat*varargin{1}; % Direct interaction
    V = V*1e-7;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ V, R ] = resist(K,I,~,L,W,H,IACS)
%RESIST: calculates the PEEC matrix vector product [R]*I for a
%collection of blocks
%   The method simply sets up an internal conductivity and calculates the
%   reistance by volume.
if nargin == 6
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IACS =  1;%0.9267;
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end
SIGMA = IACS * 5.8001e7; % International Annealed Copper Standard 100% value
L = sqrt(sum(L.^2,2));
A = sqrt(sum(W.^2,2)) .* sqrt(sum(H.^2,2));
R = L./A/SIGMA;
V = R.*I;
end
end
end
% STATIC FUNCTIONS
function vecout = blkmult(Ms,vec,a)
%BLKMATMult Block matrix structure.
%   Consider M = blkdiag(Ms.....)
% or         M = blkdiag(Ms*a(1),Ms*a(2),...)
% Then M*vec  = reshape(Ms*reshape(vec)),
%      M'*vec = reshape(Ms'*reshape(vec)).
vec = reshape(vec,size(Ms,2),[]);
if nargin > 3
    vec = bsxfun(@times,vec,a(:).');
end
vecout = Ms * vec;
vecout = vecout(:);
end
% Generate the Loop matrix. Each segment must have the same number of
% filaments
function M_ = loopMat(Nsegs_,Nfils_)
    M2 = [ones(Nfils_-1,1), spdiags(-ones(Nfils_-1,1),0,Nfils_-1,Nfils_-1)];
    M3 = cell(1,sum(Nsegs_));
    for i_ = 1:sum(Nsegs_)
        M3{i_} = M2;
    end
    M_ = blkdiag(M3{:});
end
% Terminal matrix. Each segment must have the same number of
% filaments
function T_ = termMat(Nsegs_,Nfils_)
    T_ = cell(1,numel(Nsegs_));
    for i_ = 1:numel(Nsegs_)
        tmp = sparse(1,...
            1+Nfils_*(0:Nsegs_(i_)-1),1,1,Nsegs_(i_)*Nfils_);
        T_{i_} = tmp;
    end
    T_ = blkdiag(T_{:});
end
% Generate the mesh matrix for impedance extraction.
% Nsegs = number of segments. Nfils = number of filaments per segment.
function M_ = meshMat(Nsegs_,Nfils_)
    M1 = sparse(1,(0:Nsegs_-1)*Nfils_+1,1,1,Nfils_*Nsegs_);
    M2 = [ones(Nfils_-1,1), spdiags(-ones(Nfils_-1,1),0,Nfils_-1,Nfils_-1)];
    M3 = cell(1,Nfils_);
    for i_ = 1:Nsegs_
        M3{i_} = M2;
    end
    M_ = [M1; blkdiag(M3{:})];
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
% 3. Neither the name of the copyright holder nor the names of its 
% contributors may be used to endorse or promote products derived from 
% this software without specific prior written permission.
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