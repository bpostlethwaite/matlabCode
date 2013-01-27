load('Ellipse.mat'); 
close all
topoR = makerefmat('RasterSize', size(magdata), ...
   'Latlim', [-90 90], 'Lonlim', [0 360]);

% Set up Robinson proj
figure; axesm('eqaconicstd','Origin',[lat0 lon0])

% Specify a 10x20 cell graticule
spacing = [180 360];

% Display data mapped to the graticule
h = meshm(magdata,topoR,spacing);

% Set DEM color map
demcmap('inc',magdata,10)
hold on
plotm(reflat,reflon)

plotm(reflat2,reflon2)

% demcmap(color,Z,spec) uses the color string to define a colormap. If the string is set to 'size',
% spec is the length of the colormap. If it is set to 'inc', spec is the size of the altitude range assigned to each color.
% If omitted, color is 'size' by default.