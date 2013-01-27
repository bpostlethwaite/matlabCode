function apply

global mag2data

%%
% apply gives the chosen susceptibility value to chosen prisms.
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
sus_celle=mag2data.model.sus_celle;
ind_cella=mag2data.model.ind_cella;
sus=mag2data.model.sus;
CX=mag2data.model.CX;
CZ=mag2data.model.CZ;
x_centrocella=mag2data.model.x_centrocella;
z_centrocella=mag2data.model.z_centrocella;
xspace=mag2data.model.xspace;
zspace=mag2data.model.zspace;

mag2data.controls.model.immagine=1;  % if the mask already exists, I
% mustn't do it again!

%% Application of values
%%
% sus is a vector which length is CTOT (= number of pixels = unknown 
% parametres)
sus(ind_cella)=sus_celle;
% ---------- Plot the mask ---------- %
matrix_celle=reshape(sus,CX,CZ);
imagesc(x_centrocella(1,:),z_centrocella(:,1),matrix_celle'), shading flat, grid on
set(gca, 'XAxisLocation','top', 'XTick',xspace, 'YTick',zspace, 'Position',[0.05 0.1 0.87 0.44]);
h=colorbar;
set(h, 'Position',[0.93 0.1 0.03 0.45]);
xlabel('Distance (m)','VerticalAlignment','baseline')
ylabel('Depth (m)','VerticalAlignment','middle')
axis('ij')
% ---------- Delete variable ---------- %
mag2data.model.ind_cella=[];

mag2data.model.sus=sus;
mag2data.model.matrix_celle=matrix_celle;

set(findobj('Tag','selectprismsmodelbutton'),'Enable','on')
set(findobj('Tag','selectareamodelbutton'),'Enable','on')
set(findobj('Tag','deselectmodelbutton'),'Enable','on')
set(findobj('Tag','savemodelbutton'),'Enable','on')
set(findobj('Tag','okmodelbutton'),'Enable','on')