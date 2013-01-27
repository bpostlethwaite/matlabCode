%%
% MAG2DATA is a program for 2D (21/2D) forward modelling and 2D (21/2D) 
% inversion of magnetic data. There is also a tool for simple and fast data
% filtering.
%
% Reference (1): Last, B. J. J. and Kubik, K., 1983. Compact gravity
% inversion. Geophysics, 48, 713-721.
% Reference (2): Menke, W., 1984. Geophysical Data Analysis: Discrete
% Inverse Theory. Academic Press Inc., San Diego, California, 285pp.
% Reference (3): Telford, W. M., Gedart, L. P. and Sheriff, R. E., 1990.
% Applied Geophysics. Cambridge University Press, p. 95.
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


ccc;

global mag2data % global variable

% set(0,'DefaultAxesFontSize',12);

mag2data.controls.kern='F2D';
mag2data.controls.Risoluzione=0;
mag2data.controls.angoloBeta=[];
mag2data.controls.griglia=0;
mag2data.controls.pesi=0;
mag2data.controls.initialmodel=0;
mag2data.controls.annulla=0;

menu_main;