% This file is part of the fastlitz package by Richard Y. Zhang.
% http://web.mit.edu/ryz/www
function showGeom(geom,ax)
%SHOWGEOM Render the 3D geometry modeled in geom.
if nargin == 1
    clf;
    ax = gca;
end
for ii = 1:numel(geom.O)
    showSegs(geom.O{ii},geom.L{ii},geom.W{ii},geom.H{ii},ax);
end
axis(ax,'equal','tight');
view(ax,3);

set(gcf,'Renderer','zbuffer')
material shiny
view(63,16)
camlight('left');
end

function showSegs(O,L,W,H,ax)
%SHOWSEGS
%       7----8
%      /|   /|
%     / 5----6
%    / /  / /
%   / /  / /
%  3-/--4 /
%  |/   |/
%  1----2
%
% Make eight vertices for each filament, find the convex hull for the
% entire segment of many filaments.
SIDE = [1,0.6,0.6]; CONTACT = [1, 200/255, 200/255];

if nargin == 1 && numel(O) == 4 %if input is a cell
    L = O{2};
    W = O{3};
    H = O{4};
    O = O{1};
end
% let the input be a cell array of matrices, each with 3 columns
Nsegs = numel(L);
X = cell(1,Nsegs); Y = cell(1,Nsegs); Z = cell(1,Nsegs);
X2 = cell(1,2); Y2 = cell(1,2); Z2 = cell(1,2);
for ii = 1:Nsegs
    Nfils = size(O{ii},1);
    % Compute all vertices
    vert = [O{ii};O{ii}+W{ii};O{ii}+H{ii}+W{ii};...
        O{ii}+H{ii};... Front four points
        O{ii}+L{ii};O{ii}+L{ii}+W{ii};O{ii}+L{ii}+H{ii};...
        O{ii}+L{ii}+W{ii}+H{ii};]; % Rear four points
    K = convhulln(vert).';
    
    % Divide K into surface triangles and contact trangles
    front = all(K < 4*Nfils); % First 4 Nfils vertices are all at the front
    rear = all(K > 4*Nfils); % Final 4 Nfils vertices are all at the front
    sides = ~(front|rear);
    
    % Select the patches corresponding to the sides
    x = vert(:,1); y = vert(:,2); z = vert(:,3);
    X{ii} = x(K(:,sides));
    Y{ii} = y(K(:,sides));
    Z{ii} = z(K(:,sides));
    
    % Select the patches corresponding to the front and rear
    if ii == 1
        X2{1} = x(K(:,front));
        Y2{1} = y(K(:,front));
        Z2{1} = z(K(:,front));
    end
    if ii == Nsegs
        X2{2} = x(K(:,rear));
        Y2{2} = y(K(:,rear));
        Z2{2} = z(K(:,rear));
    end
end

% Display the computed patches
patch([X{:}],[Y{:}],[Z{:}],SIDE,'linestyle','none','parent',ax);
patch([X2{:}],[Y2{:}],[Z2{:}],CONTACT,'linestyle','none','parent',ax);
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