function select_area

global mag2data

%%
% select_area selects an area of prisms.
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


set(findobj('Tag','selectprismsmodelbutton'),'Enable','off')
set(findobj('Tag','modelmenu'),'CurrentAxes',findobj('Tag','assi_griglia'));

if mag2data.controls.model.immagine==0       % Creo la maschera se non esiste
    mag2data.model.matrix_celle=reshape(mag2data.model.sus,mag2data.model.CX,mag2data.model.CZ);
    image(mag2data.model.x_centrocella(1,:),mag2data.model.z_centrocella(:,1),mag2data.model.matrix_celle'); 
    shading flat, grid on
    set(gca, 'XAxisLocation','top', 'XTick',mag2data.model.xspace, 'YTick',mag2data.model.zspace)
    xlabel('Distance (m)','VerticalAlignment','baseline')
    ylabel('Depth (m)','VerticalAlignment','middle')
    axis('ij')
end

[mag2data.model.area_in,xi,yi]=roipoly;

mag2data.model.ind_cella=find(mag2data.model.area_in'==1);

set(findobj('Tag','applymodelbutton'),'Enable','on')