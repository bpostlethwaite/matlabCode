function saveXF

global mag2data

%%
% saveXF saves processed (filtered) data in a .txt file.
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


%% Header
head1='X[m]';
head2='Data[nT or nT/m]';
header=strcat(head1, '\t', head2);

%% Saving file
[FileName,PathName] = uiputfile('*.txt','Save MAG File');
local_path=pwd;
if PathName~=0
    cd(PathName);
    fid = fopen(FileName,'w');
    fprintf(fid, header, 'newline', 'pc');
    fclose(fid);
    dlmwrite(FileName, mag2data.fileout, 'newline','pc', 'delimiter','\t', '-append', 'roffset',1, 'precision','%.3f');
    cd(local_path);
end