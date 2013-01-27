close all
clear all

load -mat campbell.mat

% Preliminaries:  Plot the grids

figure(1)
axis([min(lon) max(lon) min(lat) max(lat)])
imagesc(lon,lat,ojtopo)
set(gca,'ydir','normal')
xlabel('Longitude')
ylabel('Latitude')
title('Topography in meters')
shading interp
view(0,90)
colorbar

figure(2)
axis([min(lon) max(lon) min(lat) max(lat)])
imagesc(lon,lat,ojgrav)
set(gca,'ydir','normal')
xlabel('Longitude')
ylabel('Latitude')
title('(Free Air) Gravity in mgal')
shading interp
view(0,90)
colorbar


figure(3)
axis([min(lon) max(lon) min(lat) max(lat)])
imagesc(lon,lat,ojgeoid)
set(gca,'ydir','normal')
xlabel('Longitude')
ylabel('Latitude')
title('(Reference Geoid anomaly) Gravity in mgal')
shading interp
view(0,90)
colorbar
