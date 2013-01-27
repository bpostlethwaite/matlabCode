clear all
close all
N = 100;
h = 1;

for ii = 1:N
    for jj = 1:N
        Z(ii,jj) = h * abs(  (min([ii,jj,abs(N-ii),abs(N-jj)]) - 1)   );
    end
end

x = linspace(-1,1,100);
[X,Y] = meshgrid(x,x);

% surf(X,Y,Z)
% shading interp
% light('Position',[.5 .2 4])
% colormap copper
% 
% phi = 20*(X.^2 - Y.^4) + 50;
% 
% 
% figure(2)
% surf(X,Y,phi)
% shading interp
% light('Position',[.5 .2 .8])
% colormap winter

w = 0.402527148826027;
N = 10;
h = 1/(N+1);
G = numgrid('S',N+2);
D = delsq(G) - diag(w^2*ones(N^2,1));
[x,y] = meshgrid(1:N,1:N);
lam = 4 - 2*(cos(x*pi*h) + cos(y*pi*h)) - w^2;
lam = sort(lam(:));
E = eig(D);

[lam,E]