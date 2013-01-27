function saveinv

global mag2data

%%
% saveinv saves the results of the inversion, both profile and model.
%
% Written by Stefano Stocco on December, 2007

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

%---------- Result without the "mask signal"-----------
%% Header
head1='X[m]';
head2='Data[nT or nT/m]';
header=strcat(head1, '\t', head2);

%% Saving file
[FileName,PathName] = uiputfile('*.txt','Save calculated data (result without using the "mask signal")');
local_path=pwd;
if PathName~=0
    cd(PathName);
    fid = fopen(FileName,'w');
    fprintf(fid, header, 'newline', 'pc');
    fclose(fid);
    dlmwrite(FileName, [mag2data.Xi mag2data.Fcalcolato], 'newline','pc', 'delimiter','\t', '-append', 'roffset',1, 'precision',6);
    cd(local_path);
end

%% Header
head1='X[m]';
head2='Z[m]';
head3='Susceptibility[SI]';
header=strcat(head1, '\t', head2, '\t', head3);
 
%% Saving file
[FileName,PathName] = uiputfile('*.pam','Save model (result without using the "mask signal")');
local_path=pwd;
if PathName~=0
    cd(PathName);
    fid = fopen(FileName,'w');
    fprintf(fid, header, 'newline', 'pc');
    fclose(fid);
    dlmwrite(FileName, [mag2data.coor_centrocella mag2data.S], 'newline','pc', 'delimiter','\t', '-append', 'roffset',1, 'precision',12);
    cd(local_path);
end

%---------- Result with the "mask signal"-----------
%% Header
head1='X[m]';
head2='Data[nT or nT/m]';
header=strcat(head1, '\t', head2);

%% Saving file
[FileName,PathName] = uiputfile('*.txt','Save calculated data (result using the "mask signal")');
local_path=pwd;
if PathName~=0
    cd(PathName);
    fid = fopen(FileName,'w');
    fprintf(fid, header, 'newline', 'pc');
    fclose(fid);
    dlmwrite(FileName, [mag2data.Xi mag2data.Fcalcolatouti], 'newline','pc', 'delimiter','\t', '-append', 'roffset',1, 'precision',6);
    cd(local_path);
end

%% Header
head1='X[m]';
head2='Z[m]';
head3='Susceptibility[SI]';
header=strcat(head1, '\t', head2, '\t', head3);

%% Saving file
[FileName,PathName] = uiputfile('*.pam','Save model (result using the "mask signal")');
local_path=pwd;
if PathName~=0
    cd(PathName);
    fid = fopen(FileName,'w');
    fprintf(fid, header, 'newline', 'pc');
    fclose(fid);
    dlmwrite(FileName, [mag2data.coor_centrocella mag2data.Suti], 'newline','pc', 'delimiter','\t', '-append', 'roffset',1, 'precision',12);
    cd(local_path);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------- All the results within an ASCII file-----------

header1=strcat('X[m]', '\t', 'Z[m]', '\t', 'Susceptibility[SI] ("entire signal")', '\t', 'Susceptibility[SI] ("useful signal")');
header2=strcat('X[m]', '\t', 'Experimental Data[nT or nT/m]', '\t', 'Computed Data[nT or nT/m] ("entire signal")', '\t', 'Computed Data[nT or nT/m] ("useful signal")');

%% Saving file
[FileName,PathName] = uiputfile('*.res','Save inversion results');
local_path=pwd;
if PathName~=0
    cd(PathName);
    fid = fopen(FileName,'w');
    fprintf(fid,'INVERSION RESULT\n\n\n%%%% Inversion parameters %%%%\n\n');
    fprintf(fid,'Max iteration = %i\n',mag2data.Max_iter);
    fprintf(fid,'Noise/Signal = %i\n',mag2data.l_0);
    fprintf(fid,'Susceptibility contrast = %i\n\n',mag2data.target);
    
    fprintf(fid,'\n%%%% Inversion result ("entire signal") %%%%\n\n');
    fprintf(fid,'Perturbation number = %i\n',mag2data.eps_sus);
    fprintf(fid,'Iteration = %i\n',mag2data.min_k);
    fprintf(fid,'L2norm = %i\n\n',mag2data.L2norm);
    
    fprintf(fid,'\n%%%% Inversion result ("useful signal") %%%%\n\n');
    fprintf(fid,'Perturbation number = %i\n',mag2data.eps_susuti);
    fprintf(fid,'Iteration = %i\n',mag2data.min_kuti);
    fprintf(fid,'L2norm = %i\n\n\n\n',mag2data.L2normuti);
    
    fprintf(fid,header1);
    fclose(fid);
    dlmwrite(FileName, [mag2data.coor_centrocella mag2data.S mag2data.Suti], 'delimiter','\t', '-append', 'roffset', 1, 'precision',12);
    
    fid = fopen(FileName,'a');
    fprintf(fid,'\n\n\n');
    fprintf(fid,header2);
    fclose(fid);
    dlmwrite(FileName, [mag2data.Xi mag2data.F_obs_int mag2data.Fcalcolato mag2data.Fcalcolatouti], 'delimiter','\t', '-append', 'roffset', 1, 'precision',6);

    cd(local_path);
end