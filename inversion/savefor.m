function savefor

global mag2data

%%
% savefor saves a file with the coordinates and the susceptibility of the
% prisms. The model can be used as the initial model for an inversion.
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

%---------- Result -----------
%% Header
head1='X[m]';
head2='Z[m]';
head3='Susceptibility[SI]';
header=strcat(head1, '\t', head2, '\t', head3);
 
%% Saving file
[FileName,PathName] = uiputfile('*.pam','Save model');
local_path=pwd;
if PathName~=0
    cd(PathName);
    fid = fopen(FileName,'w');
    fprintf(fid, header, 'newline', 'pc');
    fclose(fid);
    dlmwrite(FileName, [mag2data.model.coor_centrocella mag2data.model.sus], 'newline','pc', 'delimiter','\t', '-append', 'roffset',1, 'precision',12);
    cd(local_path);
end