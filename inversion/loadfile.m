function loadfile

global mag2data

%% Load Files.
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


mag2data.controls.annulla=0;

[FileName,PathName] = uigetfile('*.txt','Select file');
LocalPath=pwd;
if PathName==0
    mag2data.controls.annulla=1;
else    
    cd(PathName);
    dati_mag=dlmread(FileName, '\t', 1, 0);
    cd(LocalPath);

    xx=dati_mag(:,1);
    mag2data.Fobs=dati_mag(:,2); 
    mag2data.lp=max(xx)-min(xx);
    mag2data.N_dati_ori=length(xx);
    mag2data.xx=xx;
    mag2data.dati_mag=dati_mag;
    
    set(findobj('Tag','modelbutton'),'Enable','on')
    set(findobj('Tag','inv2dresbutton'),'Enable','on');
    set(findobj('Tag','inv2dokbutton'),'Enable','on');
end