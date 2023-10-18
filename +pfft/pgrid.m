% This file is part of the matlab-pfft package by Richard Y. Zhang.
% http://web.mit.edu/ryz/www
classdef pgrid < handle
    %GRID grid(d,N,O)
    % To do: Add documentation at the top
    %   ALL INDEXING START AT 1. [1, 1, 1] is the origin!!!
    
    properties
        d = [10 10 10] * 1e-3; % meters
        N = 2.^[3 3 3]; 
        dims = 3;
        orig = [0 0 0];
        cnr = -3.5*10*1e-3*[1 1 1]; % bottom corner where index is (0,0,0)
        
        % Allocation padding. 
        pad = 0;
        ubnd; lbnd;
    end
    
    methods 
%==========================================================================
% CONSTRUCTOR
% Make a grid to fit points within the bounds provided, allowing for some
% padding.
%==========================================================================
function h = pgrid(delta,bnds,pad)
        if nargin == 0, return; end
        if nargin == 2, pad =0; end
        if nargin >3
            error('grid: Only 3 parameters in constructor');
        end

        % Error checks
        h.dims = size(bnds,2);
        assert(size(bnds,1)==2, 'bnds must have exactly two rows, with one upper bound and one lower bound');
        assert(numel(delta)==1 || numel(delta)==h.dims,...
            'delta must have either 1 element or dims elements');
        assert(numel(pad)==1,...
            'pad must be a scalar');
        
        % Compute the origin and span
        h.ubnd = max(bnds,[],1); h.lbnd = min(bnds,[],1);
        span = h.ubnd - h.lbnd; h.orig = h.lbnd + 0.5*span;
        
        % It can be shown that the following logic will always create a
        % grid that will have ubnd and lbnd as interior points
        h.N = ceil(span./delta)+1; 
        
        % Add the padding
        h.N = h.N + 2*pad;
        
        % Form "d" vector
        if numel(delta) == 1
            h.d = repmat(delta,1,h.dims);
        else
            h.d = delta;
        end
        
        h.cnr = h.orig - (h.d .* (h.N-1))*0.5;
        h.pad = pad;
    end
        
%==========================================================================
% PFFT-SPECIFIC METHODS
%==========================================================================
    % Gives the coordinate box that can be used to generate the Green's
    % function convolution matrix
    function cout = getCoordBox(hObj,rotdims)
        this_N = hObj.N; % Get bounds
        this_d = hObj.d;
        if nargin == 1
            rotdims = zeros(1,hObj.dims); % default to not rotating anything.
        end
        assert(numel(rotdims) == hObj.dims,'rotdims must have the same number of dimensions as the box');
        
        % If we rotate anything then we will both with finding the corners
        if any(rotdims)
            gbot = hObj.coord(1)./this_d;
        end
        
        % guide is the box for this dimension
        guide = cell(hObj.dims,1);
        for ii = 1:hObj.dims
            if rotdims(ii)
                % Hankel style
                temp = (2*gbot(ii)) + (0:(2*this_N(ii)-1));
                guide{ii} = this_d(ii)*temp([this_N(ii):end, 1:(this_N(ii)-1)]).';
            else
                % Toeplitz style
                guide{ii} = this_d(ii)*[0:(this_N(ii)), (this_N(ii)-1):-1:1].';
            end
        end

        % Repeat and expand each box guide to a total box
        if hObj.dims>1
            M = cell(1,hObj.dims);
            [M{:}] = ndgrid(guide{:});
            cout = M;
        else
            cout = guide;
        end
    end
        
    % Allocate a list of points to grid indices
    % Bounds check is now streamlined.
    function [g_ind, g_coord] = allocate(h, pnts)
        % Check dimensions & bounds
        assert(size(pnts,2) == h.dims,'Num of cols in points list must be the same as the dim of the grid');
        assert(all(all(bsxfun(@ge, pnts, h.lbnd))), 'All points must be within the box specified for this grid');
        assert(all(all(bsxfun(@le, pnts, h.ubnd))), 'All points must be within the box specified for this grid');

        % Allocate a grid point to each centroid
        g_ind = h.coord2ind(pnts);

        if nargout == 2
            % Give a list of stencil coordinates
            g_coord = h.coord(g_ind);
        end
    end
    
    function [Psparse] = proj(g,pnt,type)
    %PROJ Performs point-wise projection of a set of points onto a lattice
    %grid using Langrangian Interpolants.
    % Inputs:
    %       g   - a pgrid object used to define the grid being projected onto
    %       pnt - points to be projected
    %       type - [optional] set to true if you want a fast matrix of
    %       ones, to compute num neighbors
    % Output:
    %       Sparse form of the projection matrix.
    % Contains checks to reroute 
        if nargin == 2
            type = false;
        end
        if isstruct(pnt)
            Psparse = g.proj_bas(pnt,type);
            return
        end
        assert(isnumeric(pnt), 'Wrong input format! Either pass points or pass a basis function');
        assert(size(pnt,2) == g.dims);

        % Allocate the points to vertices or boxes. 
        % Offset the points to obtain relative distances
        [g_ind, g_coord] = g.allocate(pnt);
        pnt = bsxfun(@minus,pnt,g_coord);

        % Set up the stencil
        % Order of projectin is twice the padding.
        sten = -g.pad:g.pad;

        % Set up the first dimension
        % For each element of g_ind, displace according to sten.
        if type % just ones
            Pmat = ones(numel(sten),size(pnt,1));
        else % Lagrangian interpolants
            Pmat = lg1d(sten*g.d(1),pnt(:,1));
        end
        Pi = bsxfun(@plus,sten(:),g_ind(:).');

        % Process each subsequent dimension
        kN = [1, cumprod(g.N)]; % index stride size for each dimension
        for dim = 2:g.dims
            % Interpolate just this dimension
            if type % just ones
                tempPmat = ones(numel(sten),size(pnt,1));
            else % Lagrangian interpolants
                tempPmat = lg1d(sten*g.d(dim),pnt(:,dim));
            end
            tmpdisp = sten(:) * kN(dim); % index displacement

            % Product rule this with the previous weights
            Pmat = kronop(@times,Pmat,tempPmat);
            Pi = kronop(@plus,Pi,tmpdisp);
        end

        % Turn into a sparse matrix
        Pj = repmat(1:size(Pi,2),size(Pi,1),1);
        Psparse = sparse(Pi,Pj,Pmat,kN(end),size(Pi,2));
    end
    
    function [Psparse] = proj_bas(g,basis,type)
    %PROJ Performs point-wise projection of the quadrature points of a set
    %of bases. Basis is a struct, but MUST be formated according to a
    %specific format
    % Inputs:
    %       g   - a pgrid object used to define the grid being projected onto
    %       basis - a struct array with M elements, where M is the number
    %               of quadrature points.
    %       basis(i).weight - the scalar weight value for the i-th
    %               quadrature point
    %       basis(i).points - a N x D matrix, for the i-th quadrature point
    %               of N basis functions in D-dimensions.
    % Output:
    %       Sparse form of the projection matrix.
    if ~isstruct(basis) || ~isfield(basis,'points') || ~isfield(basis,'weight')
        error('Wrong basis input format!');
    end
        Psparse = sparse(0);
        for wi = 1:numel(basis)
            % Perform quadrature with the projection points
            Psparse = Psparse + ... 
                basis(wi).weight * ...
                g.proj(basis(wi).points,type);
        end
    end
    
    % Count the number of neighboring interactions, based on the padding
    % structure defined in this object
    function numnei = nei(g,testing, basis)
        I = g.proj(testing,true);
        if nargin == 2 || isempty(basis)
            P = I;
        else
            P = g.proj(basis,true);
        end
        
        % Non-zero pattern
        S = I' * logical(P);
        numnei = nnz(S);
    end
    
    function num = numel(g)
        num = prod(g.N);
    end
    
%==========================================================================
% CONVERSION METHODS
% Format for all conversions: one new row for each entry.
%==========================================================================
% SUBSCRIPTS <-> INDEX
    % Converts absolute subscripts to indices within this grid
    function ind = sub2ind(hObj,cin)
        % Reverse the ordering of the subscripts and turn into a cell
        cin = num2cell(cin,1);
        ind = sub2ind(hObj.N,cin{:});
    end

    % convert absolute index into subscripts within this grid
    function ss = ind2sub(hObj,ndx)
        ss = cell(1,hObj.dims);
        [ss{:}] = ind2sub(hObj.N, ndx);
        ss = cat(2,ss{:});
    end
   
%--------------------------------------------------------------------------
% COORDS -> INDEX / SUBSCRIPTS
    % convert coordinates to nearest subscript
    function cout = coord2sub(hObj,cin)
        % Assign and round
        cout = bsxfun(@minus,cin,hObj.cnr);
        cout = round(bsxfun(@rdivide,cout,hObj.d)) +1;
    end

    % convert coordinates to nearest index
    function ind = coord2ind(hObj,cin)
        ss = hObj.coord2sub(cin);
        ind = hObj.sub2ind(ss);
    end
    
%--------------------------------------------------------------------------
% INDEX / SUBSCRIPTS -> COORDS
    % convert index or subscript into absolute coordinates
    function cout = coord(hObj,cin)
        if size(cin,2) == 1 % indices
            % Convert indices into subscripts
            sub = hObj.ind2sub(cin);
        else
            sub = cin;
        end
        cout = bsxfun(@times,sub-1, hObj.d);
        cout = bsxfun(@plus,cout, hObj.cnr);
    end

    end
end

function wp = kronop(op,w1,w2)
%KRONOP Row-wise Kronecker Operator Expansion
% Applies the kronecker operator expansion along dim 1 (per-row)
% If size of dim 2 (per-column) is 1 for either inputs, will apply bsxfun
% in order to resolve the operator
%
% example 1: 
% >>w1 = [1 2;2,3]; w2 = [3,4;5,6];
% w1 = 1 2   w2 = 3 4
%      2 3        5 6
%
% >>wp = kronop(@times,w1,w2)
%    wp = 1x3 , 2x4 = 3  8
%         2x3 , 3x4   5  12
%         1x5 , 2x6   6  12
%         2x5 , 3x6   10 18
%
% example 2:
% >>w1 = [1 2;2,3]; w2 = [3;5];
% w1 = 1 2   w2 = 3
%      2 3        5
%
% >>wp = kronop(@plus,w1,w2)
%    wp = 1+3 , 2+3 = 4  5
%         2+3 , 3+3   5  6
%         1+5 , 2+5   6  7
%         2+5 , 3+5   7  8
    m1 = size(w1,1);
    m2 = size(w2,1);
    
    % Perform kronecker product expansion
    w1 = repmat(w1,m2,1);
    w2 = repmat(w2(:).',m1,1);
    w2 = reshape(w2,m1*m2,[]);
    
    % Apply operator
    wp = bsxfun(op,w1,w2);
end

function w = lg1d(x,u)
% LG1D 1-D Lagrange interpolants
% Consider the interpolant: f(uj) = \sum_i wij * f(xi)
% With x being a vector of points where f is known, and u being a vector of
% points where f must be interpolated

x = x(:); u = u(:); % Make both columns
w = 1;
for m = 1:length(x)
    hoz = u.' - x(m); % numerator
    vert = x - x(m); % denominator
    tmp = bsxfun(@rdivide,hoz,vert); % divide
    tmp(m,:) = 1; % Get rid of the m=j row.
    w = w .* tmp;
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