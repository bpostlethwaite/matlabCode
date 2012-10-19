function zeropadprofile

global mag2data

%%
% zeropedprofile zero-peds data, in order to avoid border effects.
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


% ---------- Input ---------- %
zeropadinprof=str2num(get(findobj('Tag','edit_fil_zeropadinprof'),'string'));
zeropadfineprof=str2num(get(findobj('Tag','edit_fil_zeropadfineprof'),'string'));
x=mag2data.x;
dx=mag2data.lp/(length(mag2data.Fin)-1);
% ---------- Output ---------- %
mag2data.Fout=[zeros(zeropadinprof,1); mag2data.Fin; zeros(zeropadfineprof,1)];
xoutin=[min(x)-dx:-dx:min(x)-zeropadinprof*dx]';
xoutin=flipud(xoutin);
xoutfine=[max(x)+dx:dx:max(x)+zeropadfineprof*dx]';
mag2data.xout=[xoutin; x; xoutfine];
mag2data.fileout=[mag2data.xout mag2data.Fout];
%% Saving file
saveXF;