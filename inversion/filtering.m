function filtering

global mag2data

%%
% filtering is a menu that crates variables by values setted in the grafic
% interface and shows the proper "script way" according to user's choice.
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


%% Filter type
switch mag2data.controls.filtertype
    case 'movingavarage'
        mag2data.movingavaragepoints=str2num(get(findobj('Tag','movingavaragepoints'),'string'));
        movingavaragefilter;    % Performs a moving average filter
    case 'butterworth'
        mag2data.Wp=str2num(get(findobj('Tag','edit_fil_Wp'),'string'));
        mag2data.Rp=str2num(get(findobj('Tag','edit_fil_Rp'),'string'));
        mag2data.Ws=str2num(get(findobj('Tag','edit_fil_Ws'),'string'));
        mag2data.Rs=str2num(get(findobj('Tag','edit_fil_Rs'),'string'));
        butterworthfilter;  % Performs a Butterworth filter
end
%% Spectra of the filtered signals
spectrafilt;

% ---------- Delating previous filtered data in the plot ---------- %
delete(findobj('Tag','ampfig'));
delete(findobj('Tag','dBfig'));
delete(findobj('Tag','profilofig'));
%% Plotting of the filtered data
plot_spectrafilt;