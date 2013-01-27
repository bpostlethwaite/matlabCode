function whichplot

global mag2data

%%
% whichplot switches among different plot.
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


switch mag2data.controls.computation
    case 'computing_forward'
        plot_for;
    case 'computing_inversion'
        plot_inv;
    case 'computing_synthetic'
        plot_syn;    
    case 'computing_invofsyn'
        plot_invofsyn;
        plot_iterations;
end