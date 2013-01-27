clear all
close all
%close all
n = 30;
x = -2*pi:pi/n:2*pi;
[X,Y] = meshgrid(x,x);
t = 1:0.01:10;

kx = 2*pi*X;
ky = 2*pi*Y;


for ii = 1:length(t)


Wave = cos(kx.*X./50  - t(ii)) + cos(ky.*Y./50  - t(ii));

figure(1)
surf(X,Y,Wave)
shading interp
end


