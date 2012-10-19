function grigliassi

global mag2data

%%
% grigliassi creates axes and grid for model mask.
%
% Written by Stefano Stocco on September, 2007

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
lp=mag2data.model.lp;
CX=mag2data.model.CX;
CZ=mag2data.model.CZ;
dz=mag2data.model.dz;
dx=lp/CX;
CTOT=CX*CZ;
prof=CZ*dz;
xspace=[0:dx:lp-dx/2];
zspace=[0:dz:prof-dz/2];

%% Axes and grid
set(0,'CurrentFigure',findobj('Tag','modelmenu'));
axis([0 lp 0 prof]);
set(gca,'Tag','assi_griglia')
axis('ij')
[x_centrocella,z_centrocella]=meshgrid(xspace+dx/2,zspace+dz/2);
set(gca, 'XAxisLocation','top', 'XTick',xspace, 'YTick',zspace, 'Position',[0.05 0.1 0.87 0.44])
xlabel('Distance (m)','VerticalAlignment','baseline')
ylabel('Depth (m)','VerticalAlignment','middle')
grid on
%%
% Creating two vectors with coordinate (x_centrocella,z_centrocella) into
% matrix coor_centrocella(CTOT,2)
coor_centrocella=zeros(CTOT,2);
for rig=1:CZ
    for col=1:CX
        coor_centrocella((rig-1)*CX+col,1)=x_centrocella(rig,col);
        coor_centrocella((rig-1)*CX+col,2)=z_centrocella(rig,col);
    end
end

%% Output
mag2data.model.dx=dx;
mag2data.model.CTOT=CTOT;
mag2data.model.xspace=xspace;
mag2data.model.zspace=zspace;
mag2data.model.x_centrocella=x_centrocella;
mag2data.model.z_centrocella=z_centrocella;
mag2data.model.coor_centrocella=coor_centrocella;