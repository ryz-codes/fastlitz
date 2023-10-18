% This file is part of the fastlitz package by Richard Y. Zhang.
% http://web.mit.edu/ryz/www
function [twistpath] = twist(path, mask, pitch)
%twist Twists a path with the mask as the cross section

if nargin == 0
    mask = litz.helixmask(1,5); % simple mask
    path = [linspace(0,10,100);zeros(2,100)]'; % straight path
    pitch = 10;
end

% Gather global data
nstrand = size(mask,1);

if iscell(path)
    npaths = numel(path);
    twistpath = cell(npaths,1);
    for ij = 1:npaths
        twistpath{ij} = dolitz(path{ij});
    end
    twistpath = vertcat(twistpath{:});
else
    twistpath = dolitz(path);
end

if nargout == 0
    litz.showPath(twistpath);
end

    function litzpath = dolitz(path)
        
        
        % Compute the x, y, z unit vectors for each segment of the path
        lx = path(2:end,:)-path(1:end-1,:); % Length vector
        lmag = sqrt(sum(lx.^2,2));
        ux = bsxfun(@rdivide,lx,lmag); ux = [ux;ux(end,:)];
        uy = ux*[0 1 0; -1 0 0; 0 0 0]; 
        uz = cross(ux,uy,2); 

        % Compute the parametric distance for each path segment
        % Form rotation components
        t = [0;cumsum(lmag)]/pitch;
        cost = cos(2*pi*t);
        sint = sin(2*pi*t);
        
        % Make Rx, Ry according to the following
        % [Rx] = [cos -sin] [x]
        % [Ry] = [sin  cos] [y]
        % Let Rx, Ry be matrices of nsegs x nstrands
        Rx = bsxfun(@times,cost,mask(:,1)') ...
            - bsxfun(@times,sint,mask(:,2)');
        Ry = bsxfun(@times,sint,mask(:,1)') ...
            + bsxfun(@times,cost,mask(:,2)');

        % Compute each litz path
        litzpath = cell(nstrand,1);
        for ii = 1:nstrand
            litzpath{ii} = path + bsxfun(@times,uy,Rx(:,ii)) ...
                + bsxfun(@times,uz,Ry(:,ii));
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
