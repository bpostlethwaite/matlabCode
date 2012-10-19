function butterworthfilter

global mag2data

%%
% butterworthfilter performs, of course, a butterworth filter with the
% parameters specified by the user, using the buttord, butter and filtfilt 
% Matlab function.
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


Wp=mag2data.Wp./mag2data.kny;
Ws=mag2data.Ws./mag2data.kny;
Rp=mag2data.Rp;
Rs=mag2data.Rs;
Fin=mag2data.Fin;

switch length(Wp)
    case 1
        if Wp<Ws
            ftype = 'low';
        else % Wp>Ws
            ftype = 'high';
        end
    case 2
        if Wp(1) < Ws(1) < Ws(2) < Wp(2)
            ftype = 'stop';
        else % Ws(1) < Wp(1) < Wp(2) < Ws(2)
            ftype = 'bandpass';
        end
end

% ---------- Filter ---------- %
[n,Wn]=buttord(Wp,Ws,Rp,Rs);
[b,a]=butter(n,Wn,ftype);
mag2data.Fout=filtfilt(b,a,Fin);