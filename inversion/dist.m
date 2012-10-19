function dist

global mag2data

%%
% dist marks the chosen prisms.
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


set(findobj('Tag','selectareamodelbutton'),'Enable','off')

%% Input
x_click=mag2data.model.x_click;
z_click=mag2data.model.z_click;
coor_centrocella=mag2data.model.coor_centrocella;

mag2data.controls.model.immagine=1;  % if the mask already exists, I
% mustn't do it again!

%%
% Compare coordinates (x_click,z_click) with the coordinates of the
% centre of the pixels (coor_centrocella)
for click=1:length(x_click)
    r=(((x_click(click)-coor_centrocella(:,1)).^2)+((z_click(click)-coor_centrocella(:,2)).^2)).^0.5;
    ind_cella(click)=find(r==min(r));
end
% ---------- Delete variable ---------- %
mag2data.model.x_click=[];
mag2data.model.z_click=[];

%% Output
mag2data.model.ind_cella=ind_cella;

set(findobj('Tag','applymodelbutton'),'Enable','on')