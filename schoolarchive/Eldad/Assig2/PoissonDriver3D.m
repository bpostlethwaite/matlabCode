%3D Poisson Driver
clear all
close all

%% Set constants

nn = 3:20;
L = 1;

for ii = 1:length(nn)
n1 = nn(ii); % Number of cells in x1 direction
n2 = 4*nn(ii); % Number of cells in x2 direction
n3 = nn(ii); % Number of cells in x3 direction

h1 = L/n1; % cell length in x1 direction
h2 = L/n2;
h3 = L/n3;

Lx1 = n1*h1;
Lx2 = n2*h2;
Lx3 = n3*h3;

% For now just using h and n (square cells and domain), not h1 h2 h3 so:
%h = h1;


%% Set up cells, mesh, Operators %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MESH
[x,y,z] = ndgrid(0:h1:Lx1,0:h2:Lx2,0:h3:Lx3); % Cell nodes
[xc,yc,zc] = ndgrid(h1/2:h1:Lx1-h1/2, h2/2:h2:Lx2-h2/2, h3/2:h3:Lx3-h3/2); % Cell centres
[xdx,ydx,zdx] = ndgrid(0:h1:Lx1, h2/2:h2:Lx2-h2/2, h3/2:h3:Lx3-h3/2); %Staggered in x1 cell wall x2 x3
[xdy,ydy,zdy] = ndgrid(h1/2:h1:Lx1-h1/2, 0:h2:Lx2, h3/2:h3:Lx3-h3/2); %Staggered in x2 cell wall x1 x3
[xdz,ydz,zdz] = ndgrid(h1/2:h1:Lx1-h1/2, h2/2:h2:Lx2-h2/2, 0:h3:Lx3); %Staggered in x3 cell wall x1 x2


%% Fields
%
a = 1e-6; b = pi/2 + 0.01;
u = @(t)(atan(a*t - 1/2*a).*t - 1/2*atan(a*t - 1/2*a) ...
    -1/(2*a)*log((a*t - 1/2*a).^2 + 1) + b*t);
up = @(t)(atan(a*(t-0.5)) + b);
m = @(t)(1./up(t).*exp(t));
J = @(t)(up(t).*m(t));
q = @(t)(exp(t));
%}

%{
a = 1e-6; b = pi/2 + 0.01;
u = @(t)(atan(a*t - 1/2*a).*t - 1/2*atan(a*t - 1/2*a) ...
    -1/(2*a)*log((a*t - 1/2*a).^2 + 1) + b*t);
up = @(t)(atan(a*(t-0.5)) + b);
m = @(x,y,z)(1./(up(x) + up(y) + up(z)) .*exp(x + y + z));
J = @(x,y,z)((up(x) + up(y) + up(z)).*m(t));
q = @(x,y,z)(exp(x + y + z));
%}

%{
    f = sin(pi*xc).*sin(pi*yc).*sin(pi*zc);
    fpp = -3*pi^2*sin(pi*xc).*sin(pi*yc).*sin(pi*zc);
%}

%%%%%%%%%%%%%%%%%% OPERATORS
%%%%%%%
% derivatives on walls (going from centre to nodes) <GRAD>
ddxn = @(m,k) 1/k*spdiags([-ones(m+1,1) ones(m+1,1)],[-1,0],m+1,m); 
% derivative function (calculating from nodes to centre) <DIV>
ddxc = @(m,k) 1/k*spdiags([-ones(m+1,1) ones(m+1,1)],[0,1],m,m+1); 
% Average function (Go from centres to walls (will also need ghost))
av = @(m) 0.5*spdiags([ones(m,1) ones(m,1)],[-1,0],m+1,m); 
%%%%%%% 

I1   = speye(n1);     % Create Identities of Appropriate size
I2   = speye(n2);     % Create Identities of Appropriate size
I3   = speye(n3);     % Create Identities of Appropriate size

Dn1  = ddxn(n1,h1);  % Create 1D Operators 
Dn2  = ddxn(n2,h2);  % Create 1D Operators 
Dn3  = ddxn(n3,h3);  % Create 1D Operators 

Dc1  = ddxc(n1,h1);  % Create 1D Operators 
Dc2  = ddxc(n2,h2);  % Create 1D Operators 
Dc3  = ddxc(n3,h3);  % Create 1D Operators 

Av1  = av(n1);
Av2  = av(n2);
Av3  = av(n3);


Av1([1,end]) = 0*[1,1];
Av2([1,end]) = 0*[1,1];
Av3([1,end]) = 0*[1,1];

%
Dn1([1,end]) = [2/h1; -2/h1];
Dn2([1,end]) = 0*[2/h2; -2/h2]; % Hardwire in the 0 gradient
Dn3([1,end]) = 0*[2/h3; -2/h3];
%}

%{
Dn1([1,end]) = [2/h1; -2/h1];
Dn2([1,end]) = [2/h2; -2/h2]; % Hardwire in the 0 gradient
Dn3([1,end]) = [2/h3; -2/h3];
%}

% u Boundary Conditions
%
bx = sparse(n1+1,1);
bx(1)  = 2/h1 * u(0);
bx(end) = -2/h1 * u(1);

by = sparse(n2+1,1);   % Add no boundary in for the y and z direction
bz = sparse(n3+1,1);
%}

%Bx = sparse(n1+1,n2,n3);
%{
Boundx = zeros(n1+1,n2,n3);
Boundy = zeros(n1,n2+1,n3);
Boundz = zeros(n1,n2,n3+1);

Boundx(1,:,:)   = 2/h1*(u(0) + (u(ydx(1,:,:)) + u(zdx(1,:,:))));
Boundx(end,:,:) = -2/h1*(u(1) + (u(ydx(end,:,:)) + u(zdx(end,:,:))));
Boundy(:,1,:)   = 2/h1*(u(0) + (u(xdy(:,1,:)) + u(zdy(:,1,:))));
Boundy(:,end,:) = -2/h1*(u(1) + (u(xdy(:,end,:)) + u(zdy(:,end,:))));
Boundz(:,:,1)   = 2/h1*(u(0) + (u(xdz(:,:,1)) + u(ydz(:,:,1))));
Boundz(:,:,end) = -2/h1*(u(1) + (u(xdz(:,:,end)) + u(ydz(:,:,end))));
%}

% If we change the function to be dependent on y and z will need B2 and B3
% but right now we can just reuse them since we are skipping u(y) and u(z)
% M average boundary conditions 
% m(anistropy constant, up(direction we are differentiating),x,y,z)
% 
%     
%     c1 = 1; c2 = 2;bx c3 = 3; % Anisotropy factors
    mbx = sparse(n1+1,1);   % Create sticks of boundary conditions for m
    mby = sparse(n2+1,1); 
    mbz = sparse(n3+1,1); 
    mbx(1)   = 0.5 * (1/m(-h1/2) + 1/m(h1/2));
    mbx(end) = 0.5 * (1/m(1+h1/2) + 1/m(1-h1/2));
    mby(1)   = 0.5 * (1/m(-h2/2) + 1/m(h2/2));
    mby(end) = 0.5 * (1/m(1+h2/2) + 1/m(1-h2/2));
    mbz(1)   = 0.5 * (1/m(-h3/2) + 1/m(h3/2));
    mbz(end) = 0.5 * (1/m(1+h3/2) + 1/m(1-h3/2));
% 

%%% Ramp up to 2D then to 3D
%% 3D
DN1 = kron(I3,kron(I2,Dn1));
DN2 = kron(I3,kron(Dn2,I1));
DN3 = kron(Dn3,kron(I2,I1));

DC1 = kron(I3,kron(I2,Dc1));
DC2 = kron(I3,kron(Dc2,I1));
DC3 = kron(Dc3,kron(I2,I1));

A1 = kron(I3,kron(I2,Av1));
A2 = kron(I3,kron(Av2,I1));
A3 = kron(Av3,kron(I2,I1));

Bx = kron(ones(n3,1),kron(ones(n2,1),bx));
By = kron(ones(n3,1),kron(by,ones(n1,1)));
Bz = kron(bz,kron(ones(n2,1),ones(n1,1)));

mBx = kron(ones(n3,1),kron(ones(n2,1),mbx));
mBy = kron(ones(n3,1),kron(mby,ones(n1,1)));
mBz = kron(mbz,kron(ones(n2,1),ones(n1,1)));

a1 = A1*(1./m(xc(:))) + mBx;
a2 = A2*(1./m(yc(:))) + mBy;
a3 = A3*(1./m(zc(:))) + mBz;

Am1   = spdiags(a1,0,size(A1,1),size(A1,1));
Am2   = spdiags(a2,0,size(A2,1),size(A2,1));
Am3   = spdiags(a3,0,size(A3,1),size(A3,1));
Ainv1 = spdiags(1./a1,0,size(A1,1),size(A1,1));
Ainv2 = spdiags(1./a2,0,size(A2,1),size(A2,1));
Ainv3 = spdiags(1./a3,0,size(A3,1),size(A3,1));


% Create Large operators
DIV  = [DC1,DC2,DC3];
GRAD = [DN1;DN2;DN3];
B    = [Bx;By;Bz];
S    = blkdiag(Am1,Am2,Am3);
Sinv = blkdiag(Ainv1,Ainv2,Ainv3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Solve

    Jyzero = zeros(numel(J(ydy)),1);  
    Jzzero = zeros(numel(J(zdz)),1);  
    Jvec = [J(xdx(:)) ; Jyzero(:) ; Jzzero(:) ];
 
    A = DIV*Sinv*GRAD;
    b = q(xc(:)) + DIV*Sinv*B;
    uN = A\b;
    
    
    % Residuals
    r1 = (GRAD*u(xc(:)) - B - S*Jvec);
    r2 = DIV*Jvec - q(xc(:));
    r3 = u(xc(:)) - uN;
    
    r1norm(ii) = norm(r1,'inf'); %#ok<*SAGROW>
    r2norm(ii) = norm(r2,'inf');
    r3norm(ii) = norm(r3,'inf');
    
     
end

%
figure(1)
    semilogy(nn,r1norm,'r', nn,r2norm,'b', nn,r3norm,'g', nn,1./(nn),':', nn, 1./(nn).^2,':')
    legend('r1 {grad*u - mJ}','r2 {divJ - q}','r3 {uN - uA}','O(h)','O(h^2)')
    title('Convergence Residual showing O(h^2)')
    ylabel('log residual')
    xlabel('N')
 %}

%{
    res = reshape(uN,n1,n2,n3);
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
%}









