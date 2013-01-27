clear all;
close all;

load l90_200.br;
tmp=reshape(l90_200,361,181);
tmp=flipud(tmp');
lon=0:1:360;
lat=-90:1:90;

figure(1)
imagesc(lon,lat,tmp); set (gca,'ydir','normal'); title('Langlais l90/h200'); colorbar

[Plg,Plt]=meshgrid(lon,lat);            %need grids of lat lon to make the map
figure(2);
m_proj('Azimuthal Equal-Area','longitude',209,'latitude',66,'rad',180); %Azimuthal equal-area
m_pcolor(Plg,Plt,tmp);                 %syntax for plotting map
shading flat;colormap(jet);             %colormap and shading
h=colorbar('EastOutside');
set(get(h,'title'),'string','Langlais at 200 km');    %label the colorbar
set(gca,'FontSize',12);                             %font size for labels
