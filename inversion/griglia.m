function griglia

global mag2data

%%
% griglia creates a grid with "meshgrid" function for plotting data with
% smart axes
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

mag2data.controls.griglia=1;

mag2data.dx=mag2data.lp/mag2data.CX;
mag2data.CTOT=mag2data.CX*mag2data.CZ;

mag2data.prof=mag2data.CZ*mag2data.dz;
mag2data.xspace=[0:mag2data.dx:mag2data.lp-mag2data.dx/2];
mag2data.zspace=[0:mag2data.dz:mag2data.prof-mag2data.dz/2];

[mag2data.x_centrocella,mag2data.z_centrocella]=meshgrid(mag2data.xspace+mag2data.dx/2,mag2data.zspace+mag2data.dz/2);

% Creating two vectors with coordinates (x_centrocella, z_centrocella) into
% matrix coor_centrocella(CTOT,2)
mag2data.coor_centrocella=[];                           
mag2data.coor_centrocella=zeros(mag2data.CTOT,2);

for rig=1:mag2data.CZ
    for col=1:mag2data.CX
        mag2data.coor_centrocella((rig-1)*mag2data.CX+col,1)=mag2data.x_centrocella(rig,col);
        mag2data.coor_centrocella((rig-1)*mag2data.CX+col,2)=mag2data.z_centrocella(rig,col);
    end
end