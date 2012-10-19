function lunghezza_profilo

global mag2data

%%
% lunghezza profilo writes the value "length of profile" into the inversion
% menu.
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

if mag2data.controls.annulla==0
    
    uicontrol(...
       'Parent',findobj('Tag','inversionmenu'),...
       'Tag','edit_inv_lp',...
       'Visible','on',...
       'Style','edit',...
       'Units','normalized',...
       'Position',[.275 .64 .1 .055],...
       'Enable','off',...
       'String',sprintf('%.0f',mag2data.lp),...
       'Backgroundcolor',[1 1 1]);
end
