function resmodsyn

global mag2data

%%
% resmodsyn is the core for computing resolution in the synthetic data
% process. It defines both variables and program path.
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


mag2data.controls.Risoluzione=0;
mag2data.controls.computation='computing_synthetic';

mag2data.N_dati_ori=str2num(get(findobj('Tag','edit_syn_N_dati'),'string'));
mag2data.lp=str2num(get(findobj('Tag','edit_syn_lp'),'string'));
mag2data.order=str2num(get(findobj('Tag','edit_syn_order'),'string'));
mag2data.CX=str2num(get(findobj('Tag','edit_syn_CX'),'string'));
mag2data.CZ=str2num(get(findobj('Tag','edit_syn_CZ'),'string'));
mag2data.dz=str2num(get(findobj('Tag','edit_syn_dz'),'string'));
mag2data.F=str2num(get(findobj('Tag','edit_syn_F'),'string'));
mag2data.I=str2num(get(findobj('Tag','edit_syn_I'),'string'));
mag2data.Beta=str2num(get(findobj('Tag','edit_syn_Beta'),'string'));
mag2data.hsi=str2num(get(findobj('Tag','edit_syn_hsi'),'string'));
mag2data.hss=str2num(get(findobj('Tag','edit_syn_hss'),'string'));
mag2data.Qn=str2num(get(findobj('Tag','edit_syn_Qn'),'string'));
if mag2data.Qn~=0
    mag2data.I_R=str2num(get(findobj('Tag','edit_syn_I_R'),'string'));
    mag2data.Beta_R=str2num(get(findobj('Tag','edit_syn_Beta_R'),'string'));
else
    mag2data.I_R=0;
    mag2data.Beta_R=0;
end
mag2data.controls.campo_tot=1;
set(findobj('Tag','syntheticmenu'),'Visible','off');

whichkernel;
resolution;
griglia;         % crea una griglia per plottare i dati con assi "furbi"
plot_resolution;     % plotto res e R
disp('Press any key to continue!!');
pause
set(findobj('Tag','syntheticmenu'),'Visible','on');   % riapro la maschera dell'inversione per procedere
mag2data.controls.Risoluzione=1;