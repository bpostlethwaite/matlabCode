function plot_spectra

global mag2data

%%
% plot_spectra plots spectra evaluated in spectra function. Y axes are both
% in amplitude and dB; X axes are in wave number and wave length. It also
% plots the signal.
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


%% Evaluation of wave lengths from wave numbers
passo=round(max(mag2data.k))/10;
ktick=[min(mag2data.k):passo:max(mag2data.k)];
lamtick=(2*pi)./ktick;
lamtick=round(lamtick*10)/10;
clc

%% Plots
% ---------- Amplitude spectrum ---------- %
set(findobj('Tag','spectrummenu'),'CurrentAxes',findobj('Tag','assi_amp'));
plot(mag2data.k,mag2data.Samp), grid on
ylabel('Amplitude')
set(gca,'XLim',[min(mag2data.k) max(mag2data.k)], 'XTick',ktick)
set(gca,'Tag','assi_amp')
hold on

% ---------- dB spectrum ---------- %
set(findobj('Tag','spectrummenu'),'CurrentAxes',findobj('Tag','assi_dB'));
plot(mag2data.k,mag2data.SdB), grid on
xlabel('Wave number k (rad/m)')
ylabel('dB')
set(gca,'XLim',[min(mag2data.k) max(mag2data.k)], 'XTick',ktick, 'XAxisLocation','top')
set(gca,'Tag','assi_dB')
hold on

% ---------- Wave lengths axes ---------- %
set(findobj('Tag','spectrummenu'),'CurrentAxes',findobj('Tag','assi_lam'));
plot(mag2data.k,zeros(1,length(mag2data.k)),'black')
xlabel('Wave length \lambda (m)')
set(gca,'XLim',[min(mag2data.k) max(mag2data.k)], 'XTick',ktick, 'XTickLabel', lamtick)
set(gca,'Tag','assi_lam')

% ---------- Signal ---------- %
set(findobj('Tag','spectrummenu'),'CurrentAxes',findobj('Tag','assi_profilo'));
plot(mag2data.x,mag2data.Fin), grid on
xlabel('Distance (m)')
ylabel('Magnetic anomaly (nT - nT/m)')
set(gca,'XLim',[min(mag2data.x) max(mag2data.x)])
set(gca,'Tag','assi_profilo')
hold on