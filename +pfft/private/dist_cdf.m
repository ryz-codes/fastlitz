% This file is part of the matlab-pfft package by Richard Y. Zhang.
% http://web.mit.edu/ryz/www
function [dist,cdf] = dist_cdf(testing, basis, numsamp)
%DIST_CDF Summary of this function goes here
%   Detailed explanation goes here

% ALL NORMS ARE L_INF

% Obtain centroids and radiuses
[cent,radt,numt] = get_dat(testing);
if ~isempty(basis)
    [cenb,radb,numb] = get_dat(basis);
else
    cenb = cent; radb = radt; numb = numt;
end

% Set up interaction pairs
numpair = numb*numt;
numsamp = min(numpair,numsamp);
k = randi(numpair,numsamp,1); % April 2014 - Sample with replacement
[i,j] = ind2sub([numt,numb],k);

% Compute the distances
d = cent(i,:) - cenb(j,:);
dist = max(abs(d),[],2) - radt(i) - radb(j); % Triangle inequality
dist = max(dist,0); % Floor at 0. Distance cannot be negative.

% Sort and compute cdf
dist = sort(dist);
cdf = (1:numsamp)*(numpair / numsamp);
cdf = reshape(cdf,size(dist));

% Pull out the unique values
[dist,ia] = unique(dist);
cdf = cdf(ia);
end

function [cen,rad,num,dims] = get_dat(dat)
if isnumeric(dat)
    cen = dat;
    rad = 0;
    num = size(dat,1);
    dims = size(dat,2);
elseif isstruct(dat) && isfield(dat,'points') && isfield(dat,'weight')
    % Merge all the data points 
    pnts = cat(3,dat.points); % One slice for every quad point
    cen = mean(pnts,3); % Centroids for every basis
    rad = bsxfun(@minus,pnts,cen); % distance from each quad point to centroid
    rad = max(max(rad,[],2),[],3); % max(L_inf(radius))
    num = size(pnts,1);
    dims = size(pnts,2);
else
    error('Wrong input format for testing function!');
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