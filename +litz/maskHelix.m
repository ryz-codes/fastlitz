% This file is part of the fastlitz package by Richard Y. Zhang.
% http://web.mit.edu/ryz/www
function [mask] = maskHelix(srad, nstrand)
%maskHelix Generates a helix-style mask for a set of strands.
% This mask basically puts the strands in a tight circle
% Inputs:
%   srad - the strands will be spaced exactly srad apart.
%   nstrand - the number of strands

if nargin == 0 % Defaults
    srad = 1;
    nstrand = 9;
end

if nstrand == 1 % One strand takes no further processing
    mask = [0 0];
else
    % Radius of the circle made by squeezing nstrand smaller circles together.
    % We compute this radius by simple geometry arguments, on the isosceles
    % triangle whose edge lengths are hrad, hrad and 2*srad respectively.
    hrad = srad / sin(pi/nstrand);

    % We rotate around the circle using a parametric argument
    t = linspace(0,2*pi,nstrand+1);
    t = t(1:end-1);
    mask = hrad*[cos(t);sin(t)]';
end

if nargout == 0
    litz.showMask(mask,srad);
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