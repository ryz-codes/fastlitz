% This file is part of the fastlitz package by Richard Y. Zhang.
% http://web.mit.edu/ryz/www
function showPEEC(geom,ax)
%SHOWGEOM Show the 3D geometry modeled in geom.
if nargin == 1
    ax = gca;
end
O2 = expand(geom.O);
L2 = expand(geom.L);
W2 = expand(geom.W);
H2 = expand(geom.H);
showFils(O2,L2,W2,H2,ax);
end

function output = expand(input)
output = vertcat(input{:});
output = vertcat(output{:});
end

function [ output_args ] = showFils(O,L,W,H,ax)
%SHOWFILS Summary of this function goes here
%   Detailed explanation goes here
%       7----8
%      /|   /|
%     / 5----6
%    / /  / /
%   / /  / /
%  3-/--4 /
%  |/   |/
%  1----2
%
% Make eight vertices
%
% If the number of fils to draw exceeds 200, we won't draw the front and 
% back faces.

draw_face = size(O,1) < 200;
draw_face = true;
if nargin == 1 && numel(O) == 4 %if input is a cell
    L = O{2};
    W = O{3};
    H = O{4};
    O = O{1};
end

vert{1} = O;
vert{2} = O+W;
vert{3} = O+H;
vert{4} = O+W+H;
vert{5} = O+L;
vert{6} = O+L+W;
vert{7} = O+L+H;
vert{8} = O+L+W+H;

% Make six faces
if draw_face
    [f1x f1y f1z] = getface (1,2,4,3);
    [f6x f6y f6z] = getface (5,6,8,7);
end
[f2x f2y f2z] = getface (1,5,7,3);
[f3x f3y f3z] = getface (2,6,8,4);
[f4x f4y f4z] = getface (3,4,8,7);
[f5x f5y f5z] = getface (1,2,6,5);


% Draw the faces
cla(ax);
for ii = 1:size(O,1)
    if draw_face
        patch(f1x(ii,:), f1y(ii,:), f1z(ii,:), [1 0.7 0.7], 'parent', ax); %pink
        patch(f6x(ii,:), f6y(ii,:), f6z(ii,:), [1 0.7 0.7], 'parent', ax); %pink
    end
    
    patch(f2x(ii,:), f2y(ii,:), f2z(ii,:), 'red', 'parent', ax);
    patch(f3x(ii,:), f3y(ii,:), f3z(ii,:), 'red', 'parent', ax);
    patch(f4x(ii,:), f4y(ii,:), f4z(ii,:), 'red', 'parent', ax);
    patch(f5x(ii,:), f5y(ii,:), f5z(ii,:), 'red', 'parent', ax);
    
%     % Uncomment below to have more realistic pics
%     patch(f2x(ii,:), f2y(ii,:), f2z(ii,:), 'red','linestyle','none');
%     patch(f3x(ii,:), f3y(ii,:), f3z(ii,:), 'red','linestyle','none');
%     patch(f4x(ii,:), f4y(ii,:), f4z(ii,:), 'red','linestyle','none');
%     patch(f5x(ii,:), f5y(ii,:), f5z(ii,:), 'red','linestyle','none');
end
view(ax,3);
axis(ax,'equal');
% light('Position',[1 1 1]);
% material dull
% lighting phong
%shading interp

    % Private function to sort out the vertices
    function [fx, fy, fz] = getface (a,b,c,d)
        fx = [vert{a}(:,1) vert{b}(:,1) vert{c}(:,1) vert{d}(:,1)];
        fy = [vert{a}(:,2) vert{b}(:,2) vert{c}(:,2) vert{d}(:,2)];
        fz = [vert{a}(:,3) vert{b}(:,3) vert{c}(:,3) vert{d}(:,3)];
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
