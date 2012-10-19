function showprofile

global mag2data

%%
% showprofile opens a file with two or three vectors (x, total magnetic 
% field, vertical gradient) and plots them.
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


loadfile;
pause(0.1)

if mag2data.controls.annulla==0

    dati_mag=mag2data.dati_mag;

    [r,c]=size(dati_mag);

    x=dati_mag(:,1);
    Ftot=dati_mag(:,2);

    if c==2
        figure
        plot(x,Ftot)
        legend('Total field (bottom sensor) or Vertical gradient')
        xlabel('Distance (m)')
        ylabel('Magnetic anomaly (nT - nT/m)')
        set(gca,'XLim',[0 max(x)])
    else    % c=3
        Gradvert=dati_mag(:,3);
        figure
        plot(x,Ftot,'b',x,Gradvert,'m')
        legend('Total field (bottom sensor)','Vertical gradient')
        xlabel('Distance (m)')
        ylabel('Magnetic anomaly (nT - nT/m)')
        set(gca,'XLim',[0 max(x)])
    end
end