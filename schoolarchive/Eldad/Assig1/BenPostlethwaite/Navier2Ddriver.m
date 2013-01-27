% Marker and Cell - Navier.
clear all
close all

test = 1;

L = 1;
nn = 4:4:96;

for ii = 1:length(nn)
n = nn(ii);
h = L/n;

[xc,yc]   = ndgrid(h/2:h:L-h/2, h/2:h:L-h/2); % Cell centres
[xdx,ydx] = ndgrid(h:h:L-h, h/2:h:L-h/2); % x walls 
[xdy,ydy] = ndgrid(h/2:h:L-h/2, h:h:L-h); % y walls 

%% Operators
% 2nd order double derivative centred difference operator
ddxx = @(m,k) 1/k^2*spdiags([ones(m+1,1) -2*ones(m+1,1) ones(m+1,1)],[-1,0,1],m-1,m-1); 

% derivative function (calculating from centres to walls) 
ddxc = @(m,k) 1/k*spdiags([-ones(m+1,1) ones(m+1,1)],[0,1],m-1,m); 

% derivative function (calculating from walls to centres)
ddxw = @(m,k) 1/k*spdiags([-ones(m+1,1) ones(m+1,1)],[-1,0],m,m-1); 

% Averager for plotting vectors
av = @(m) 0.5*spdiags([ones(m,1) ones(m,1)],[-1,0],m,m-1);

% Boundary conditions on the average matrix
Av = av(n);
Av([1,end]) = 1;

%% Scale to 2D
I = speye(n,n);
Dxx = kron(I,ddxx(n,h));
Dyy = kron(ddxx(n,h),I);
Dpx = kron(I,ddxc(n,h));
Dpy = kron(ddxc(n,h),I);
Dx = kron(I,ddxw(n,h));
Dy = kron(ddxw(n,h),I);
S1 = sparse(size(Dxx,2),size(Dxx,1)); % sparse zero matrix for super matrix
S2 = sparse(size(Dpx,2),size(Dpx,2));
Avx = kron(I,Av);
Avy = kron(Av,I);

%f1 = 9.8*xdx;
%f2 = zeros(n-1,n);
%f2(1:end) = 10;

%% Define Force Fields, 
%Here is pure curl (get spiral flow)
f1 = 0.5 - ydx;
f2 =  - 0.5 + xdy;

grav = 3;
f1 = f1 + grav;
% Centre force fields for plotting
fc1 = reshape(Avx*f1(:),size(xc));
fc2 = reshape(Avy*f2(:),size(yc));

%% TESTING
a = 0.01;
ux    = (1 - 2*ydx) .* exp( -( ((xdx - 0.5).^2)/a + ((ydx - 0.5).^2)/a) );
uy    = -(1 - 2*xdy) .* exp( -( ((xdy - 0.5).^2)/a + ((ydy - 0.5).^2)/a) );
pres  = 2*(xc);
F1    = ( (2*ydx-1).*(2*a - (1 - 2*xdx).^2) .* exp( -( ((xdx - 0.5).^2)/a + ((ydx - 0.5).^2)/a) ) )./a^2 + 2;
F2    = -( (2*xdy-1).*(2*a - (1 - 2*ydy).^2) .* exp( -( ((xdy - 0.5).^2)/a + ((ydy - 0.5).^2)/a) ) )./a^2;

%% Matrix solns and reformat

Z  = zeros(n,n);

% Build Super matrix
A = [Dxx, S1, Dpx
     S1, Dyy, Dpy
     Dx, Dy,  S2 ];
 
 
if test == 1
    f1 = F1;
    f2 = F2;
end


b = [f1(:);f2(:);Z(:)];

X = A\b;

N = length(xdx(:));
u1 = reshape(X(1:N),size(xdx));
u2 = reshape(X(N+1:2*N),size(xdy));
P  = reshape(X(2*N+1 : end),size(xc));

if test == 1
    resdiff = reshape(Dx*ux(:) + Dy*uy(:),size(xc));
    
    res1x = abs(u1 - ux);
    res1y = abs(u2 - uy);
    res2  = abs(pres - P);

    r1(ii) = norm(res1x,'inf');
    r2(ii) = norm(res1y,'inf');
    r3(ii) = norm(resdiff,'inf');
end

end

% Average velocity components and centre for plotting
fc1 = reshape(Avx*f1(:),size(xc));
fc2 = reshape(Avy*f2(:),size(yc));
uc1 = reshape(Avx*u1(:),size(xc));
uc2 = reshape(Avy*u2(:),size(yc));


%% Plotting


figure(1)
semilogy(nn,r1,'--',nn,r2,'-.',nn,r3,nn,1./nn,':',nn,1./nn.^2,':')
legend('convergence ux','convergence uy','du/dx - du/dy','O(h)','O(h^2)')
title(sprintf('Convergence in velocity solutions (infinite norm)\n and derivative residual'))

figure(2)
quiver(xc,yc,fc1,fc2)
title('Input Forcing')

figure(3);
quiver(xc,yc,uc1,uc2)
title('Computed Vector Velocity Field')

figure(4)
imagesc(P)
colorbar
title('Computed Pressure gradient')

