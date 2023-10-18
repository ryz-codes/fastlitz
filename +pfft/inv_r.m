% This file is part of the matlab-pfft package by Richard Y. Zhang.
% http://web.mit.edu/ryz/www
% 
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
%INV_R vectorized 1/r green's function, with zero at the origin.
% Please input the arguments as a cell, one cell for each dimension.
% 
% Has an internal setting on the minimum distance before setting things to
% zero.
function g = inv_r(X)
MIN_DIST = eps;
Z_DIST_SQ = 0;

assert(iscell(X));

dist = Z_DIST_SQ;
for ii = 1:numel(X)
    dist = dist + X{ii}.^2;
end
dist = sqrt(dist);
g = zeros(size(dist));

% Find the non-zero points
z_idx = dist > MIN_DIST;
g(z_idx) = 1./dist(z_idx);