function deselect

global mag2data

%%
% deselect deselects the chosen prisms, giving to them the value 0.
%
% Written by Stefano Stocco on September, 2007
% Modified by Stefano Stocco on December, 2008

% Copyright (C) 2007 Stefano Stocco
% 
% This file is part of MAG2DATA.
% 
% MAG2DATA is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 3 of the License, or (at your
% option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License along
% with this program; if not, see <http://www.gnu.org/licenses>.


%% Input
x_des=mag2data.model.x_des;
z_des=mag2data.model.z_des;
coor_centrocella=mag2data.model.coor_centrocella;
x_centrocella=mag2data.model.x_centrocella;
z_centrocella=mag2data.model.z_centrocella;
matrix_celle=mag2data.model.matrix_celle;
sus=mag2data.model.sus;
CX=mag2data.model.CX;
CZ=mag2data.model.CZ;
xspace=mag2data.model.xspace;
zspace=mag2data.model.zspace;

%%
% Compare coordinates (x_des,z_des) with the coordinates of the
% centre of the pixels (coor_centrocella)
for des=1:length(x_des)
    r=(((x_des(des)-coor_centrocella(:,1)).^2)+((z_des(des)-coor_centrocella(:,2)).^2)).^0.5;
    ind_des(des)=find(r==min(r));
end
% ---------- Application of values 0 ---------- %
sus(ind_des)=0;
% ---------- Plot the mask ---------- %
matrix_celle=reshape(sus,CX,CZ);
imagesc(x_centrocella(1,:),z_centrocella(:,1),matrix_celle'), shading flat, grid on
set(gca, 'XAxisLocation','top', 'XTick',xspace, 'YTick',zspace, 'Position',[0.05 0.1 0.87 0.44])
h=colorbar;
set(h, 'Position',[0.93 0.1 0.03 0.45]);
xlabel('Distance (m)','VerticalAlignment','baseline')
ylabel('Depth (m)','VerticalAlignment','middle')
axis('ij')
% ---------- Delete variable ---------- %
mag2data.model.x_des=[];
mag2data.model.z_des=[];
mag2data.model.ind_des=[];

%% Output
mag2data.model.sus=sus;
mag2data.model.matrix_celle=matrix_celle;

set(findobj('Tag','applymodelbutton'), 'Enable','on')