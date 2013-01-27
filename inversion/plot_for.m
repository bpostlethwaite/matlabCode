function plot_for

global mag2data

%%
% plot_for plots the result of the forward modelling.
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
Xi=mag2data.Xi;
xx=mag2data.xx;
Fobs=mag2data.Fobs;
x_centrocella=mag2data.model.x_centrocella;
z_centrocella=mag2data.model.z_centrocella;
matrix_celle=mag2data.model.matrix_celle*10^3;
xspace=mag2data.model.xspace;
zspace=mag2data.model.zspace;

set(findobj('Tag','modelmenu'),'CurrentAxes',mag2data.hh);

if mag2data.controls.campo_tot==1
    plot(xx,Fobs,'.', Xi,mag2data.F_inf,'r', 'MarkerSize',8, 'LineWidth',1.5);
    legend('Observed total field','Computed total field')
    ylabel('Magnetic anomaly (nT)')
    set(gca, 'XLim',[min(Xi) max(Xi)], 'XTickLabel',[])
else
    plot(xx,Fobs,'.', Xi,mag2data.GradV,'r', 'MarkerSize',8, 'LineWidth',1.5);
    legend('Observed vertical gradient','Computed vertical gradient')
    ylabel('Magnetic anomaly (nT/m)')
    set(gca,'XLim',[min(Xi) max(Xi)], 'XTickLabel',[])
end