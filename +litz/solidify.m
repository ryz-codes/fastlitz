% This file is part of the fastlitz package by Richard Y. Zhang.
% http://web.mit.edu/ryz/www
function [geom] = solidify(path,wrad,drad,odd)
%SOLIDIFY Takes an input of path (or cell of path), and outputs a circular 
% wire that traces the path.
%                     __________
%  ---------   =>    O__________)
%
% wrad: wire radius
% drad: wire cross-sectional discretization

order = 2; % Keep this on 2 unless otherwise

if nargin ==0
    path = [0 0 0; 1 0 0];
    wrad = 1;
    drad = 2;
end
if nargin < 4
    odd = true;
end

% Construct the cross-section
% make a non-uniform grid.
%tmp = linspace(-1,1,2*drad);
if odd
    tmp = linspace(-1,1,drad*2);
else
    tmp = linspace(-1,1,drad*2+1);
end
tmp = abs(tmp).^(1/order) .* sign(tmp);
[gi, gj] = ndgrid(tmp(1:end-1));
[gw, gh] = ndgrid(diff(tmp));

% Compute the distances of the centroids
dist = (gi+gw/2).^2+(gj+gh/2).^2;

% Sort by distance from center, delete all the filaments with distances 
% greater than 1 in order to make a circle
[sdist,idx] = sort(dist(:));
idx(sdist>1) = [];
gi=gi(idx); gj=gj(idx); gw=gw(idx); gh=gh(idx);

numfils = numel(gi);
geom.mask = [gi(:),gw(:),gj(:),gh(:)]*wrad;


if iscell(path) %multiple paths passed in as cells
    numpaths = numel(path);
    geom.O = cell(numpaths,1); geom.L = cell(numpaths,1);
    geom.W = cell(numpaths,1); geom.H = cell(numpaths,1);
    for ij = 1:numpaths
        [geom.O{ij}, geom.L{ij}, ...
         geom.W{ij}, geom.H{ij}] = dopath(path{ij});
    end
else %just a single path
    [geom.O{1}, geom.L{1}, ...
     geom.W{1}, geom.H{1}] = dopath(path);
end

if nargout ==0
    litz.showPEEC(geom);
    view(120,20);
end

    function [O,L,W,H] = dopath(path)
    % Converts a single path into filaments    
        
        Nsegs = size(path,1)-1;
        
        % Compute the x, y, z unit vectors for each segment of the path
        lx = path(2:end,:)-path(1:end-1,:); % Length vector
        [~,ux] = vecnorm(lx);  % Get the lengths
        [~,uw] = vecnorm(ux*[0 1 0; -1 0 0; 0 0 0]); % Rotate on xy plane
        uh = cross(ux,uw,2);
        
        wy = uw*wrad; % Width vector
        hz = uh*wrad; % Height vector
        numsegs = size(lx,1);

        % Each segment must be packaged inside a cell
        O = cell(numsegs,1); L = cell(numsegs,1); 
        W = cell(numsegs,1); H = cell(numsegs,1);
        for ii = 1:numsegs % each segment
            tmp = gi(:)*wy(ii,:) + gj(:)*hz(ii,:);
            O{ii} = bsxfun(@plus,path(ii,:),tmp);
            L{ii} = repmat(lx(ii,:),numfils,1);
            W{ii} = gw(:)*wy(ii,:);
            H{ii} = gh(:)*hz(ii,:);
        end
    end
end

function [norms normvecs] = vecnorm(vecs)
% each row is a vector
norms = sqrt(sum(vecs.^2,2));
if nargout == 2
normvecs = bsxfun(@rdivide,vecs,norms);
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