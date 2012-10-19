function spectrafilt

global mag2data

%%
% spectrafilt evaluates the spectra of a filtered signal in wave number 
% domain using the Matlab fft function. The outputs are spectra in
% amplitude and dB.
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


%% Data in: x, signal, length of profile
x=mag2data.x;
Fout=mag2data.Fout;
lp=mag2data.lp;
dx=mag2data.lp/(length(mag2data.Fout)-1);

%% Spectrum
if round(length(Fout)/2)==length(Fout)/2
    zero_ped=length(Fout);
else
    zero_ped=length(Fout)+1;
end
kny=pi/dx;  % Nyquist wave number
deltak=2*pi/((zero_ped-1)*dx);
k=0:deltak:kny; % Wave number vector
% ---------- Fast fourier transform ---------- %
S=fft(Fout,zero_ped);
Samp=abs(S(1:length(k)));
SdB=-20*log10(max(Samp)./Samp);
% ---------- Output ---------- %
mag2data.kny=kny;
mag2data.k=k;
mag2data.Samp=Samp;
mag2data.SdB=SdB;