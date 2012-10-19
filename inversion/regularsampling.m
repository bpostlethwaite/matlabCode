function regularsampling

global mag2data

%%
% regularsampling creates evenly sampling data from unevenly sampling data
% using interp1 Matlab function and save them.
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


% ---------- Input ---------- %
dx=str2num(get(findobj('Tag','edit_fil_dx'),'string'));
x=mag2data.x;
xint=[min(x):dx:round(max(x))]';
% ---------- Output ---------- %
mag2data.Fout=interp1(x,mag2data.Fin,xint,'linear','extrap');
mag2data.fileout=[xint mag2data.Fout];

%% Saving file
saveXF;

