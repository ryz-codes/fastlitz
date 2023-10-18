% This file is part of the fastlitz package by Richard Y. Zhang.
% http://web.mit.edu/ryz/www
function [mask] = maskFilled(srad,ns)
%MASKFILLED Mask that fills a circular bundle with strands 
% Inputs:
%         srad - the radius of a single strand
%         wrad - the radius of the wire bundle
if nargin == 0
    srad = 1;
    ns = 5;
end

% Approx radius strands
nr = sqrt(ns)*pi;
% Number of layers approx sqrt(3)/2 * radius
mask = hexgrid(srad,round(2*nr));

if ns < 7
    mask = litz.maskHelix(srad,ns);
    return;
end

while size(mask,1) > ns
    rad2 = sum(mask.^2,2);
    [~,ri] = max(rad2);
    rm = mask(ri,:);
    mask(ri,:) = []; %delete
    mask = bsxfun(@plus,mask,rm/size(mask,1));
end

if nargout == 0
    litz.showMask(mask,srad);
    disp(size(mask,1));
end
end

function mask = hexgrid(d,n)
[gx,gy] = ndgrid(1:n);
% 45 degree rotation
gx2 = gx-gy;
gy2 = gx+gy;

% Remove offset
mask = [gx2(:),gy2(:)];
mask = bsxfun(@minus,mask,mean(mask));

% Scale both directions
mask = bsxfun(@times,mask,d*[1,sqrt(3)]);
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