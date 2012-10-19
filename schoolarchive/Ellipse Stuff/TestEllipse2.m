%TEST ELLIPSE II
%results
close all
clear all
%%
load l90_200.br; % Loads data
load('BenColormaps','ellipsecolor'); %Load colormap
magdata=reshape(l90_200,361,181); % Reshapes data
magdata=flipud(magdata'); % Reshape data
lon=0:1:360;
lat=-90:1:90;
% [Plg,Plt]=meshgrid(lon,lat);            %need grids of lat lon to make the map

%% Main Ellipse control
R = 1000;
semimajor = 0.5*pi*R*cosd(60);
%semiminor =  semimajor*0.2;
lat0 = 89.9; % Centre lat
lon0 = 0; %Centre Lon
offset = 0; % Degrees of Offset for ellipse
ecc = 0;
mainellipse = [semimajor,ecc];   % Defining parameters for main ellipse 
units = 'degrees'; % Some necessary parameter, might change when I find out what it does

%% Ellipsoid body on which we will plot our ellipse
mars_radius = R;
mars_eccentricity = 0;
planet = [mars_radius,mars_eccentricity];

%% Compute lon & lat of ellipse
[reflat,reflon] = ellipse1(lat0,lon0,mainellipse,offset,[],planet); %Large ellipse lat and lon

%% Distance Minimization

minorellipse = [0.5*semimajor,ecc]; %Second ellipse for computing normals
num_pts(1)=20;  % number pf points in outer ellipse
num_pts(2) = 9; % number of points in inner ellipse
az2 = [0:10:350;2:10:352]'; % Small inner angle (For track spacing)
az1 = [az2(:,1)-8,az2(:,2)+8]; % Larger Outer angle (For track spacing)

rng = -60;  % Degree range for the profile tracks
[prof_lat prof_lon] = getprofiles(lat0,lon0,mainellipse,minorellipse,offset,az1,az2,planet,units,num_pts,rng);


%% MESH AND GRID
topoR = makerefmat('RasterSize', size(magdata), ... % Make a graticule mesh
   'Latlim', [-90 90], 'Lonlim', [0 360]);  % I really don't know what this does yet
spacing = [180 360]; % Set spacing for graticule (should be around data grid size)

%% Figure 1, the square data
figure(1)
imagesc(lon,lat,magdata); set (gca,'ydir','normal'); title('Langlais l90/h200'); colorbar
colormap(ellipsecolor)

%% Figure 2, Globe projection
figure(2)

axesm ('ortho', 'Frame', 'on', 'Grid', 'on','Origin',[0 lon0],'MeridianLabel','on','MLabelParallel','equator','ParallelLabel','on','PLabelMeridian','prime');
title('Test - Normal Profile Test')
plotm(prof_lat,prof_lon)
plotm(reflat,reflon)

%% Figure 3 - Play with this figure to get best results
figure(3)
axesm ('eqdazim', 'Frame', 'on', 'Grid', 'on','Origin',[lat0 lon0],'MeridianLabel','on','MLabelParallel','equator','ParallelLabel','on','PLabelMeridian','prime');
plotm(reflat,reflon)
plotm(prof_lat,prof_lon)
title('Test - Normal Profile Test')

hold off