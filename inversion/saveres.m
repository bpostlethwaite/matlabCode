function saveres

global mag2data

%%
% saveres creates a figure with the option to save the result of the
% inversion.
%
% Written by Stefano Stocco on December, 2008

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

%% Save Result figures
figure(...
    'Tag','saveres',...
    'Menubar','none',...
    'Numbertitle','off',...
    'Units','normalized',...
    'Color',[0.5235    0.5785    0.7471],...
    'Name','Save result',...
    'Position',[0.35 0.4 0.3 0.2],...
    'RendererMode','manual');
saveres=findobj('Tag','saveres');

%% Buttons
uicontrol(...
   'Parent',saveres,...
   'Visible','on',...
   'Style','text',...
   'Units','normalized',...
   'Position',[.05 .75 .9 .1],...
   'Enable','on',...
   'String','Do you want to save the result of the inversion?',...
   'Backgroundcolor',[0.5235    0.5785    0.7471]);
uicontrol(...
   'Parent',saveres,...
   'Style','pushbutton',...
   'Units','normalized',...
   'Position',[.2 .35 .2 .2],...
   'Enable','On',...
   'String','Yes',...
   'Callback','close(findobj(''Tag'',''saveres'')); saveinv; set(findobj(''Tag'',''inversionmenu''),''Visible'',''On'')');
uicontrol(...
   'Parent',saveres,...
   'Style','pushbutton',...
   'Units','normalized',...
   'Position',[.6 .35 .2 .2],...
   'Enable','On',...
   'String','No',...
   'Callback','close(findobj(''Tag'',''saveres'')); set(findobj(''Tag'',''inversionmenu''),''Visible'',''On'')');