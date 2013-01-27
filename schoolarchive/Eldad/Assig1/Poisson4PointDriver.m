%3D Poisson Driver
clear all
close all

%% Set constants
%kk = 1;

Lx1 = 1; % Our Length Scale (Goes from 0 -> 1)
Lx2 = Lx1;
Lx3 = Lx1;
nn =  [4:20];

for ii = 1:length(nn)
n1 = nn(ii); % Number of cells in x1 direction
n2 = n1; % Number of cells in x2 direction
n3 = n1; % Number of cells in x3 direction

h1 = Lx1/n1; % cell length in x1 direction
h2 = Lx2/n2; % cell length in x2 direction
h3 = Lx3/n3; % cell length in x3 direction

% For now just using h and n (square cells and domain), not h1 h2 h3 so:
h = h1;
n = n1;

%% Build m and q

% fictitious m  (3 m's for anisotropy) THESE NEED TO BE INVERTED 1/m
% (for real use of parameters, here as 1 does not matter.
m1 = ones(n1,n2,n3);
m2 = m1;
m3 = m2;

% determine q's. This will be two delta functions one a negative one a
% positive


%% Set up cells, mesh, Operators %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MESH
[x1,x2,x3] = ndgrid(0:h1:Lx1,0:h2:Lx2,0:h3:Lx3); % Cell nodes
[xc1,xc2,xc3] = ndgrid(h1/2:h1:Lx1-h1/2,h2/2:h2:Lx2-h2/2,h3/2:h3:Lx3-h3/2); % Cell nodes
[x1f1,x2f1,x3f1] = ndgrid(h1/2:h1:Lx1-h1/2,0:h2:Lx2,0:h3:Lx3); %Staggered in x1 cell wall x2 x3
[x1f2,x2f2,x3f2] = ndgrid(0:h1:Lx1,h2/2:h2:Lx2-h2/2,0:h3:Lx3); %Staggered in x2 cell wall x1 x3
[x1f3,x2f3,x3f3] = ndgrid(0:h1:Lx1,0:h2:Lx2,h3/2:h3:Lx3-h3/2); %Staggered in x3 cell wall x1 x2

%%%%%%%%%%%%%%%%%% OPERATORS
%%%%%%%
% derivatives on walls (going from centre to nodes) <GRAD>
ddxn = @(m,k) 1/k*spdiags([-ones(m+1,1) ones(m+1,1)],[-1,0],m+1,m); 
% derivative function (calculating from nodes to centre) <DIV>
ddxc = @(m,k) 1/k*spdiags([-ones(m+1,1) ones(m+1,1)],[0,1],m,m+1); 
% Average function (Go from centres to walls (will also need ghost))
av = @(m) 0.5*spdiags([ones(m,1) ones(m,1)],[-1,0],m+1,m); 
% 4 point derivative method for going from centre to noedes <GRAD>
dd4n = @(m,k) 1/(24*k)*spdiags([ones(m+1,1),-27*ones(m+1,1),27*ones(m+1,1),-ones(m+1,1)],-2:1,m+1,m); 
% 4 point derivative method for going from centre to noedes <GRAD>
dd4c = @(m,k) 1/(24*k)*spdiags([ones(m+1,1),-27*ones(m+1,1),27*ones(m+1,1),-ones(m+1,1)],-1:2,m,m+1); 


%%%%%%% 


Ix1 = speye(n1); % Create Identities of Appropriate size
Ix2 = speye(n2);
Ix3 = speye(n3);

Dnx1  = dd4n(n1,h1);  % Create 1D Operators 
Dnx2  = dd4n(n2,h2);
Dnx3  = dd4n(n3,h3);

Dcx1  = dd4c(n1,h1);  % Create 1D Operators 
Dcx2  = dd4c(n2,h2);
Dcx3  = dd4c(n3,h3);

Avx1  = av(n1);
Avx2  = av(n2);
Avx3  = av(n3);

% correct for BC  (DIRICHLET!!!)
Dnx1(1,1:3) = [90, -20, 18/5]/(24*h1);
Dnx1(2,1:3) = [-30,28,-6/5]/(24*h1);
Dnx1(end-1,end-2:end) = [6/5,-28,30]/(24*h1);
Dnx1(end,end-2:end) = [-18/5,20, -90]/(24*h1);

Dcx1(1,1:4) = [-23, 21, 3, -1]/(24*h1);
Dcx1(end,end-3:end) = [1, -3, -21, 23]/(24*h2);


Avx1([1,end])  = 1;
Avx2([1,end])  = 1;
Avx3([1,end])  = 1;

%%% Ramp up to 2D then to 3D
%% 3D
Dn1 = kron(Ix3,kron(Ix2,Dnx1));
Dn2 = kron(Ix3,kron(Dnx1,Ix1));
Dn3 = kron(Dnx1,kron(Ix2,Ix1));

Dc1 = kron(Ix3,kron(Ix2,Dcx1));
Dc2 = kron(Ix3,kron(Dcx1,Ix1));
Dc3 = kron(Dcx1,kron(Ix2,Ix1));

A1 = kron(Ix3,kron(Ix2,Avx1));
A2 = kron(Ix3,kron(Avx2,Ix1));
A3 = kron(Avx3,kron(Ix2,Ix1));

DIV  = [Dc1,Dc2,Dc3];
GRAD = [Dn1;Dn2;Dn3];


Am1 = spdiags(A1*m1(:),0,size(A1,1),size(A1,1));
Am2 = spdiags(A2*m2(:),0,size(A2,1),size(A2,1));
Am3 = spdiags(A3*m3(:),0,size(A3,1),size(A3,1));

Sinv = blkdiag(Am1,Am2,Am3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fields
f = sin(pi*xc1).*sin(pi*xc2).*sin(pi*xc3);
%fp = 2*pi*cos(2*pi*x1f1).*sin(2*pi*x2f1).*sin(2*pi*x3f1)...
%    + 2*pi*sin(2*pi*x1f2).*cos(2*pi*x2f2).*sin(2*pi*x3f2)...
%    + 2*pi*sin(2*pi*x1f3).*sin(2*pi*x2f3).*cos(2*pi*x3f3);
%fp = 2*pi*cos(2*pi*x1).*sin(2*pi*x2).*sin(2*pi*x3)...
%    + 2*pi*sin(2*pi*x1).*cos(2*pi*x2).*sin(2*pi*x3)...
%    + 2*pi*sin(2*pi*x1).*sin(2*pi*x2).*cos(2*pi*x3);
fpp = -3*pi^2*sin(pi*xc1).*sin(pi*xc2).*sin(pi*xc3);
%% Solve
% Solving A = DIV*Sinv*GRAD in 2 steps to take advantage of backslash
% operator in matlab

SG = Sinv\GRAD;
A = DIV*SG;
u = A\fpp(:);

% Residual (computed solution minus function)
relres = abs(u - f(:))./f(:); % 0.001785
res = abs(u - f(:));
%u = reshape(u,size(q));

r(ii) = norm(res,'inf');


end

figure(1)
semilogy(nn,r,'--', nn,1./(nn).^2,':', nn, 1./(nn).^3,':', nn,1./(nn).^4,':')
legend('residual','O(h^2)','O(h^3)','O(h^4)')
title('Convergence Residual showing O(h^4)')
ylabel('log residual')
xlabel('N')

res = reshape(res,size(xc1));
     figure(2)
     imagesc(squeeze(res(:,:,round(nn(ii)/2))))
     colorbar;
     title('X Y Plane')
     
     figure(3)
     imagesc(squeeze(res(:,round(nn(ii)/2),:)))
     colorbar;
     title('X Z Plane')
     
     figure(4)
     imagesc(squeeze(res(round(nn(ii)/2),:,:)))
     colorbar;
     title('Y Z Plane')
     
     figure(5)
     imagesc(squeeze(res(:,:,end)))
     colorbar;
     title('X Y Bottom surface Plane')
     

 