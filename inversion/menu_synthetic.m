function munu_synthetic

global mag2data

%%
% menu_synthetic creates the synthetic figure menu, in order to compute the
% forward modelling of a magnetic body, ealuating both total field and
% vertical gradient, with the possibility of introducing remanent
% magnetisation into the causative body. It is then possible to perform an
% inversion procedure starting from the synthetic data.
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


%% Synthetic figure
figure(...
    'Tag','syntheticmenu',...
    'Menubar','none',...
    'Visible','on',...
    'Numbertitle','off',...
    'Units','normalized',...
    'Color',[0.5235    0.5785    0.7471],...
    'Name','SYNTHETIC DATA',...
    'Position',[0.25 0.3 0.5 0.45],...
    'RendererMode','manual');
syntheticmenu=findobj('Tag','syntheticmenu');

%% Synthetic menu
% ---------- Frame 1 ---------- %
uicontrol(...
   'Parent',syntheticmenu,...
   'Visible','on',...
   'Style','frame',...
   'Units','normalized',...
   'Position',[.03 .7 .36 .25],...
   'Enable','on',...
   'Backgroundcolor',[0.5235    0.5785    0.7471]);
% ---------- Profile ---------- %
uicontrol(...
   'Parent',syntheticmenu,...
   'Visible','on',...
   'Style','text',...
   'Units','normalized',...
   'Position',[.021 .92 .09 .05],...
   'Enable','on',...
   'String','Profile',...
   'Backgroundcolor',[0.5235    0.5785    0.7471]);
uicontrol(...
   'Parent',syntheticmenu,...
   'Visible','on',...
   'Style','text',...
   'Units','normalized',...
   'Position',[.04 .85 .22 .05],...
   'Enable','on',...
   'String','Number of Data',...
   'Backgroundcolor',[0.7    0.75    0.85]);
uicontrol(...
   'Parent',syntheticmenu,...
   'Tag','edit_syn_N_dati',...
   'Visible','on',...
   'Style','edit',...
   'Units','normalized',...
   'Position',[.275 .85 .1 .058],...
   'Enable','on',...
   'String',[101],...
   'Backgroundcolor',[1 1 1]);
uicontrol(...
   'Parent',syntheticmenu,...
   'Visible','on',...
   'Style','text',...
   'Units','normalized',...
   'Position',[.04 .75 .22 .05],...
   'Enable','on',...
   'String','Profile Length (m)',...
   'Backgroundcolor',[0.7    0.75    0.85]);
uicontrol(...
   'Parent',syntheticmenu,...
   'Tag','edit_syn_lp',...
   'Visible','on',...
   'Style','edit',...
   'Units','normalized',...
   'Position',[.275 .75 .1 .058],...
   'Enable','on',...
   'String',[100],...
   'Backgroundcolor',[1 1 1]);

% ---------- Frame 2 ---------- %
uicontrol(...
   'Parent',syntheticmenu,...
   'Visible','on',...
   'Style','frame',...
   'Units','normalized',...
   'Position',[.03 .5 .36 .15],...
   'Enable','on',...
   'Backgroundcolor',[0.5235    0.5785    0.7471]);
% ---------- Data interpolation ---------- %
uicontrol(...
   'Parent',syntheticmenu,...
   'Visible','on',...
   'Style','text',...
   'Units','normalized',...
   'Position',[.021 .62 .2 .05],...
   'Enable','on',...
   'String','Data Interpolation',...
   'Backgroundcolor',[0.5235    0.5785    0.7471]);
uicontrol(...
   'Parent',syntheticmenu,...
   'Visible','on',...
   'Style','text',...
   'Units','normalized',...
   'Position',[.04 .55 .22 .05],...
   'Enable','on',...
   'String','Order of Interpolation',...
   'Backgroundcolor',[0.7    0.75    0.85]);
uicontrol(...
   'Parent',syntheticmenu,...
   'Tag','edit_syn_order',...
   'Visible','on',...
   'Style','edit',...
   'Units','normalized',...
   'Position',[.275 .55 .1 .058],...
   'Enable','on',...
   'String',[1],...
   'Backgroundcolor',[1 1 1]);

% ---------- Frame 3 ---------- %
uicontrol(...
   'Parent',syntheticmenu,...
   'Visible','on',...
   'Style','frame',...
   'Units','normalized',...
   'Position',[.03 .02 .36 .43],...
   'Enable','on',...
   'Backgroundcolor',[0.5235    0.5785    0.7471]);
% ---------- Prisms ---------- %
uicontrol(...
   'Parent',syntheticmenu,...
   'Visible','on',...
   'Style','text',...
   'Units','normalized',...
   'Position',[.021 .42 .09 .05],...
   'Enable','on',...
   'String','Prisms',...
   'Backgroundcolor',[0.5235    0.5785    0.7471]);
uicontrol(...
   'Parent',syntheticmenu,...
   'Visible','on',...
   'Style','text',...
   'Units','normalized',...
   'Position',[.04 .35 .22 .05],...
   'Enable','on',...
   'String','Prisms along x',...
   'Backgroundcolor',[0.7    0.75    0.85]);
uicontrol(...
   'Parent',syntheticmenu,...
   'Tag','edit_syn_CX',...
   'Visible','on',...
   'Style','edit',...
   'Units','normalized',...
   'Position',[.275 .35 .1 .058],...
   'Enable','on',...
   'String',[25],...
   'Backgroundcolor',[1 1 1]);
uicontrol(...
   'Parent',syntheticmenu,...
   'Visible','on',...
   'Style','text',...
   'Units','normalized',...
   'Position',[.04 .25 .22 .05],...
   'Enable','on',...
   'String','Prisms along z',...
   'Backgroundcolor',[0.7    0.75    0.85]);
uicontrol(...
   'Parent',syntheticmenu,...
   'Tag','edit_syn_CZ',...
   'Visible','on',...
   'Style','edit',...
   'Units','normalized',...
   'Position',[.275 .25 .1 .058],...
   'Enable','on',...
   'String',[5],...
   'Backgroundcolor',[1 1 1]);
uicontrol(...
   'Parent',syntheticmenu,...
   'Visible','on',...
   'Style','text',...
   'Units','normalized',...
   'Position',[.04 .15 .22 .05],...
   'Enable','on',...
   'String','Thickness (m)',...
   'Backgroundcolor',[0.7    0.75    0.85]);
uicontrol(...
   'Parent',syntheticmenu,...
   'Tag','edit_syn_dz',...
   'Visible','on',...
   'Style','edit',...
   'Units','normalized',...
   'Position',[.275 .15 .1 .058],...
   'Enable','on',...
   'String',[2],...
   'Backgroundcolor',[1 1 1]);
uicontrol(...
   'Parent',syntheticmenu,...
   'Visible','on',...
   'Style','pushbutton',...
   'Units','normalized',...
   'Position',[.05 .05 .15 .068],...
   'Enable','on',...
   'String','Resolution',...
   'Backgroundcolor',[0.7    0.75    0.85],...
   'Callback','resmodsyn;');
uicontrol(...
    'Parent',syntheticmenu,...
    'Tag','modelbutton',...
    'Visible','on',...
    'Style','pushbutton',...
    'Units','normalized',...
    'Position',[.22 .05 .15 .065],...
    'Enable','on',...
    'String','Model',...
    'Callback','mag2data.model.lp=str2num(get(findobj(''Tag'',''edit_syn_lp''),''string'')); mag2data.model.CX=str2num(get(findobj(''Tag'',''edit_syn_CX''),''string'')); mag2data.model.CZ=str2num(get(findobj(''Tag'',''edit_syn_CZ''),''string'')); mag2data.model.dz=str2num(get(findobj(''Tag'',''edit_syn_dz''),''string'')); set(findobj(''Tag'',''syntheticmenu''),''visible'',''off''); model_syn;');

% ---------- Frame 4 ---------- %
uicontrol(...
   'Parent',syntheticmenu,...
   'Visible','on',...
   'Style','frame',...
   'Units','normalized',...
   'Position',[.61 .6 .36 .35],...
   'Enable','on',...
   'Backgroundcolor',[0.5235    0.5785    0.7471]);
% ---------- Magnetic field parameters ---------- %   
uicontrol(...
   'Parent',syntheticmenu,...
   'Visible','on',...
   'Style','text',...
   'Units','normalized',...
   'Position',[.601 .92 .13 .05],...
   'Enable','on',...
   'String','Field Data',...
   'Backgroundcolor',[0.5235    0.5785    0.7471]);
uicontrol(...
   'Parent',syntheticmenu,...
   'Visible','on',...
   'Style','text',...
   'Units','normalized',...
   'Position',[.62 .84 .22 .075],...
   'Enable','on',...
   'String','Total Field Intensity F (nT)',...
   'Backgroundcolor',[0.7    0.75    0.85]);
uicontrol(...
   'Parent',syntheticmenu,...
   'Tag','edit_syn_F',...
   'Visible','on',...
   'Style','edit',...
   'Units','normalized',...
   'Position',[.855 .85 .1 .058],...
   'Enable','on',...
   'String',[46000],...
   'Backgroundcolor',[1 1 1]);
uicontrol(...
   'Parent',syntheticmenu,...
   'Visible','on',...
   'Style','text',...
   'Units','normalized',...
   'Position',[.62 .75 .22 .05],...
   'Enable','on',...
   'String','Inclination I (�)',...
   'Backgroundcolor',[0.7    0.75    0.85]);
uicontrol(...
   'Parent',syntheticmenu,...
   'Tag','edit_syn_I',...
   'Visible','on',...
   'Style','edit',...
   'Units','normalized',...
   'Position',[.855 .75 .1 .058],...
   'Enable','on',...
   'String',[60],...
   'Backgroundcolor',[1 1 1]);
uicontrol(...
   'Parent',syntheticmenu,...
   'Visible','on',...
   'Style','text',...
   'Units','normalized',...
   'Position',[.62 .63 .22 .09],...
   'Enable','on',...
   'String','Angle Profile Direction / Magnetic North (�)',...
   'Backgroundcolor',[0.7    0.75    0.85]);
uicontrol(...
   'Parent',syntheticmenu,...
   'Tag','edit_syn_Beta',...
   'Visible','on',...
   'Style','edit',...
   'Units','normalized',...
   'Position',[.855 .65 .1 .058],...
   'Enable','on',...
   'String',[0],...
   'Backgroundcolor',[1 1 1]);

% ---------- Frame 5 ---------- %
uicontrol(...
   'Parent',syntheticmenu,...
   'Visible','on',...
   'Style','frame',...
   'Units','normalized',...
   'Position',[.61 .34 .36 .21],...
   'Enable','on',...
   'Backgroundcolor',[0.5235    0.5785    0.7471]);
% ---------- Height of Sensors ---------- %
uicontrol(...
   'Parent',syntheticmenu,...
   'Visible','on',...
   'Style','text',...
   'Units','normalized',...
   'Position',[.601 .52 .2 .05],...
   'Enable','on',...
   'String','Height of Sensors',...
   'Backgroundcolor',[0.5235    0.5785    0.7471]);
uicontrol(...
   'Parent',syntheticmenu,...
   'Visible','on',...
   'Style','text',...
   'Units','normalized',...
   'Position',[.62 .46 .22 .05],...
   'Enable','on',...
   'String','Bottom Sensor (m)',...
   'Backgroundcolor',[0.7    0.75    0.85]);
uicontrol(...
   'Parent',syntheticmenu,...
   'Tag','edit_syn_hsi',...
   'Visible','on',...
   'Style','edit',...
   'Units','normalized',...
   'Position',[.855 .46 .1 .058],...
   'Enable','on',...
   'String',[0.3],...
   'Backgroundcolor',[1 1 1]);
uicontrol(...
   'Parent',syntheticmenu,...
   'Visible','on',...
   'Style','text',...
   'Units','normalized',...
   'Position',[.62 .37 .22 .05],...
   'Enable','on',...
   'String','Top Sensor (m)',...
   'Backgroundcolor',[0.7    0.75    0.85]);
uicontrol(...
   'Parent',syntheticmenu,...
   'Tag','edit_syn_hss',...
   'Visible','on',...
   'Style','edit',...
   'Units','normalized',...
   'Position',[.855 .37 .1 .058],...
   'Enable','on',...
   'String',[1.3],...
   'Backgroundcolor',[1 1 1]);

% ---------- Frame 6 ---------- %
uicontrol(...
   'Parent',syntheticmenu,...
   'Visible','on',...
   'Style','frame',...
   'Units','normalized',...
   'Position',[.61 .02 .36 .27],...
   'Enable','on',...
   'Backgroundcolor',[0.5235    0.5785    0.7471]);
% ------------ Remanent magnetisation ------------------------
uicontrol(...
   'Parent',syntheticmenu,...
   'Visible','on',...
   'Style','text',...
   'Units','normalized',...
   'Position',[.6 .26 .27 .05],...
   'Enable','on',...
   'String','Remanent magnetisation',...
   'Backgroundcolor',[0.5235    0.5785    0.7471]);
uicontrol(...
   'Parent',syntheticmenu,...
   'Visible','on',...
   'Style','text',...
   'Units','normalized',...
   'Position',[.62 .2 .22 .05],...
   'Enable','on',...
   'String','Konigsberger ratio (%)',...
   'Backgroundcolor',[0.7    0.75    0.85]);
uicontrol(...
   'Parent',syntheticmenu,...
   'Tag','edit_syn_Qn',...
   'Visible','on',...
   'Style','edit',...
   'Units','normalized',...
   'Position',[.855 .2 .1 .058],...
   'Enable','on',...
   'String',[0],...
   'Backgroundcolor',[1 1 1]);
uicontrol(...
   'Parent',syntheticmenu,...
   'Visible','on',...
   'Style','text',...
   'Units','normalized',...
   'Position',[.62 .12 .22 .05],...
   'Enable','on',...
   'String','Vertical Angle',...
   'Backgroundcolor',[0.7    0.75    0.85]);
uicontrol(...
   'Parent',syntheticmenu,...
   'Tag','edit_syn_I_R',...
   'Visible','on',...
   'Style','edit',...
   'Units','normalized',...
   'Position',[.855 .12 .1 .058],...
   'Enable','on',...
   'String',[],...
   'Backgroundcolor',[1 1 1]);
uicontrol(...
   'Parent',syntheticmenu,...
   'Visible','on',...
   'Style','text',...
   'Units','normalized',...
   'Position',[.62 .04 .22 .05],...
   'Enable','on',...
   'String','Horizontal Angle',...
   'Backgroundcolor',[0.7    0.75    0.85]);
uicontrol(...
   'Parent',syntheticmenu,...
   'Tag','edit_syn_Beta_R',...
   'Visible','on',...
   'Style','edit',...
   'Units','normalized',...
   'Position',[.855 .04 .1 .058],...
   'Enable','on',...
   'String',[],...
   'Backgroundcolor',[1 1 1]);

% ---------- Go! Button ---------- %
uicontrol(...
   'Parent',syntheticmenu,...
   'Tag','synokbutton',...
   'Visible','on',...
   'Style','pushbutton',...
   'Units','normalized',...
   'Position',[.40 .05 .09 .08],...
   'Enable','off',...
   'String','OK',...
   'Callback','synthetic;');
uicontrol(...
   'Parent',syntheticmenu,...
   'Visible','on',...
   'Style','pushbutton',...
   'Units','normalized',...
   'Position',[.42 .85 .16 .1],...
   'Enable','on',...
   'String','Main Window',...
   'BackgroundColor',[1 0 0],...
   'Callback','close(findobj(''Tag'',''syntheticmenu'')); close(findobj(''Tag'',''modelmenu'')); set(findobj(''Tag'',''mainmenu''),''Visible'',''on'')');

% uicontrol(...
%    'Parent',syntheticmenu,...
%    'Visible','on',...
%    'Style','pushbutton',...
%    'Units','normalized',...
%    'Position',[.51 .05 .09 .08],...
%    'Enable','on',...
%    'String','Help',...
%    'Callback','set(findobj(''Tag'',''help_synfig''),''visible'',''on'')');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%