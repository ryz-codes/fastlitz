% This file is part of the fastlitz package by Richard Y. Zhang.
% http://web.mit.edu/ryz/www
function showMask(mask,srad)
    % Error check
    assert(size(mask,2) == 2, 'Mask must have 2 columns');
    
    % Plot the centers with "x"s
    plot(mask(:,1),mask(:,2),'x','MarkerSize',20);
    hdl = gca; hold on
    
    % Plot a circle for each strand
    for ii = 1:size(mask,1)
        plotcirc(hdl,mask(ii,:),srad,1);
    end
    
    % Plot the a bit circle for the bundle
    wrad = sqrt(max(sum(mask.^2,2))) + srad;
    plotcirc(hdl,[0,0],wrad,2);
    
    hold off;
    axis image
end

function plotcirc(hdl,cen,rad,wei)
    t = linspace(0,2*pi,50); % 40 sides per circle
    x = rad*cos(t)+cen(1);
    y = rad*sin(t)+cen(2);
    if nargin == 3
    plot(hdl,x,y,'k');
    else
    plot(hdl,x,y,'k','LineWidth',wei);
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