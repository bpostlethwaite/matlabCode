% A test file later to be used as a function for generating 2D function.
clear all
close all
L = 1;
n = 500;
h = L/n;
[xn,yn]   = ndgrid(0:h:L, 0:h:L); % Cell centres
[Xi,Yi]   = ndgrid(1:n+1);
Rm = randn(size(xn));

a = 1e-2;
b = [0.45,0.25,0.75,1.1];
pos = [.25, .25
       .25, .75
       .75, .25
       .75, .75];
hills = 0;
for ii = 1:length(pos)
    hills = hills + b(ii)*(exp( -( ((xn - pos(ii,1)*L).^2)./a + ((yn - pos(ii,2)*L).^2)./a  )));
end

%% hole
a = 1e-3;
pos = [0.5,0.1];
hole = -(exp( -( ((xn - pos(1)*L).^2)./a + ((yn - pos(2)*L).^2)./a  +2 )));
%% cup shape domain

b = 10;
c = 1.6;
d = 10;
cup = b*( c*(yn-L/2) ).^d;
%% saddle domain
samp = 3;
saddle = samp * ((xn - L/2).^2 - (yn - L/2).^2);
%% Plains
pos = [0.75,0.5];
rad = 0.1*(n+1);
plains = sqrt((Xi - pos(1)*(n+1)).^2 + (Yi - pos(2)*(n+1)).^2);

topo = saddle + hills + cup + hole;
topo(plains <= rad) = topo(round(pos(1)*(n+1)),round(pos(2)*(n+1)));

topo = topo + 0.01*Rm;
%% Plots
%mesh(xn,yn,bowl)


mesh(xn,yn,topo)
%shading interp
%colormap gray
