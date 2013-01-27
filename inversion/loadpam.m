function loadpam

global mag2data

%% Load Initial Model.
%
% Written by Stefano Stocco on December, 2008

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


[FileName,PathName] = uigetfile('*.pam','Select file');
LocalPath=pwd;
if PathName~=0
    cd(PathName);
    dati_mag=dlmread(FileName, '\t', 1, 0);
    cd(LocalPath);
    mag2data.coor_centrocella_initialmodel=dati_mag(:,1:2);
    mag2data.S=dati_mag(:,3);
    
    mag2data.controls.initialmodel=1;
end

