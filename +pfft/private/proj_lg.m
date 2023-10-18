% This file is part of the matlab-pfft package by Richard Y. Zhang.
% http://web.mit.edu/ryz/www
function [Pmat,Pi,Psparse] = proj_lg(g,pnt,ord)
%PROJ_LG Performs point-wise projection of a set of points onto a lattice
%grid using Langrangian Interpolants.
% Inputs:
%       g   - a pgrid object used to define the grid being projected onto
%       pnt - points to be projected
%       ord - order of projection. 
% Output:
%       Column-condensed form of the projection matrix.
if nargin < 3
    ord = 2; % Default to quadratic interpolation
end
assert(size(pnt,2) == g.dims);
assert(mod(ord,2) == 0, 'Only even orders implemented currently');

% Allocate the points to vertices or boxes. 
% Offset the points to obtain relative distances
[g_ind, g_coord] = g.allocate(pnt,ceil(ord/2));
pnt = bsxfun(@minus,pnt,g_coord);

% Set up the stencil
sten = -ord/2:ord/2;

% Set up the first dimension
% For each element of g_ind, displace according to sten.
Pmat = lg1d(sten*g.d(1),pnt(:,1));
Pi = bsxfun(@plus,sten(:),g_ind(:).');

% Process each subsequent dimension
kN = [1, cumprod(g.N)]; % index jump size for each dimension
for dim = 2:g.dims
    % Interpolate just this dimension
    tempPmat = lg1d(sten*g.d(dim),pnt(:,dim));
    tmpdisp = sten(:) * kN(dim); % index displacement
    
    % Product rule this with the previous weights
    Pmat = kronop(@times,Pmat,tempPmat);
    Pi = kronop(@plus,Pi,tmpdisp);
end

% Turn into a sparse matrix
if nargout == 3
Pj = repmat(1:size(Pi,2),size(Pi,1),1);
Psparse = sparse(Pi,Pj,Pmat,kN(end),size(Pi,2));
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