% This file is part of the fastlitz package by Richard Y. Zhang.
% http://web.mit.edu/ryz/www
function [geom] = wire(len,pitch,ns,srad,insu,drad,odd)
%LITZ.WIRE Recursive wire generation. 
% wire(len,pitch,ns,srad,insu,drad)
% len:     Length of the straight wire
% pitch:   1 x k vector of pitches, one for each level, starting from the
%          bottom
import litz.*
if nargin == 0
    len = 2e-2;
    %pitch = [0.02 -0.03 0.045];
    %ns = [5 5 5];
    pitch = [0.005, 0.1];
    ns = [5 5];
    srad = 3*1.6e-4/2;
    drad = 3;
end

if nargin < 5 || isempty(insu)
    insu = 0.08; % Default 8% insulation
end
if nargin < 6 || isempty(drad)
    drad = 2; % Default 2 cross-sectional discretization
end
if nargin < 7 || isempty(odd)
    odd = true; % Default 2 cross-sectional discretization
end

% Error checks
assert(numel(pitch) == numel(ns));
assert(numel(srad) == 1);
pitch = pitch(:); ns = ns(:);
levs = numel(pitch);

% Division per pitch
dl = min(abs([len;pitch])) / 20;
% Make the straight path
path = 0:dl:len;
path = [path; zeros(2,length(path))]';

% Compute the radius of each level
brad = zeros(levs+1,1); brad(1) = srad*(1+insu);
mask = cell(levs,1);
for ii = 1:levs
    
    thisn = ns(ii);
    if thisn>6, thisn=6;end
    
    % Compute the innermost circle radius
    wrad = fzero(@(x) getD2(x,pitch(ii),thisn)-4*brad(ii)^2,0);
    wrad = sqrt(wrad);
    
    % Compute the radius of the mask
    tmpmask = maskFilled(brad(ii),thisn);
    R = max(sqrt(sum(tmpmask.^2,2)));
    
    % Stretch factor
    fac = wrad/R;
    
    % Make mask and stretch
    mask{ii} = maskFilled(brad(ii),ns(ii))*fac;
    R = sqrt(sum(mask{ii}.^2,2));
    
    brad(ii+1) = max(R) + brad(ii);   
end

% Twist the path top down
for ii = levs:-1:1
    path = twist(path,mask{ii},pitch(ii));
end

geom = solidify(path,srad,drad,odd);

if nargout == 0
    showGeom(geom);
        view(50,20)% Check to see collision
    grid off
    axis off
end
function D2 = getD2(R2,pitch,n)
% D2 = 2R2 - 2R2cos p + P2(p-phi2)^2
% R2 sin p + P2(p-phi2) = 0

P2 = (pitch/2/pi)^2;
phi2 = 2*pi/n;

% Solve for p
if isfinite(P2)
    p = fzero(@(p) R2.*sin(p) + P2.*(p-phi2),[0,pi]);
    % Apply the equation
    D2 = 2*R2*(1-cos(p)) + P2*(p-phi2).^2;
else
    D2 = 2*R2*(1-cos(phi2));
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