function plot_spectrafilt

global mag2data

%%
% plot_spectrafilt plots spectra evaluated in spectrafilt function. Y axes 
% are both in amplitude and dB; X axes are in wave number and wave length. 
% It also plots the filtered signal.
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


%% Plots
% ---------- Amplitude spectrum ---------- %
set(findobj('Tag','spectrummenu'),'CurrentAxes',findobj('Tag','assi_amp'));
plot(mag2data.k,mag2data.Samp,'r', 'Tag','ampfig');

% ---------- dB spectrum ---------- %
set(findobj('Tag','spectrummenu'),'CurrentAxes',findobj('Tag','assi_dB'));
plot(mag2data.k,mag2data.SdB,'r', 'Tag','dBfig');

% ---------- Signal ---------- %
set(findobj('Tag','spectrummenu'),'CurrentAxes',findobj('Tag','assi_profilo'));
plot(mag2data.x,mag2data.Fout,'r', 'Tag','profilofig');