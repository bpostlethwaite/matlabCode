%3D Poisson Driver
clear all
close all

%% Set constants
%kk = 1;
%%%%%%%%%%%%%%%%
mode = 1;      %  Pick Mode 1 for working simple function. 
%%%%%%%%%%%%%%%%  Mode 2 is discontinuous function that is not working.
                  % and should allow to test derivatives.
L = 1; % Our Length Scale (Goes from 0 -> 1)

nn = 3:20;

for ii = 1:length(nn)
n = nn(ii); % Number of cells in x1 direction
%n2 = n; % Number of cells in x2 direction
%n3 = n; % Number of cells in x3 direction

h1 = L/n; % cell length in x1 direction
h2 = h1;
h3 = h1;

Lx1 = L;
Lx2 = L;
Lx3 = L;

% For now just using h and n (square cells and domain), not h1 h2 h3 so:
h = h1;


%% Set up cells, mesh, Operators %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MESH
[x,y,z] = ndgrid(0:h1:Lx1,0:h2:Lx2,0:h3:Lx3); % Cell nodes
[xc,yc,zc] = ndgrid(h1/2:h1:Lx1-h1/2, h2/2:h2:Lx2-h2/2, h3/2:h3:Lx3-h3/2); % Cell centres
[xdx,ydx,zdx] = ndgrid(0:h1:Lx1, h2/2:h2:Lx2-h2/2, h3/2:h3:Lx3-h3/2); %Staggered in x1 cell wall x2 x3
[xdy,ydy,zdy] = ndgrid(h1/2:h1:Lx1-h1/2, 0:h2:Lx2, h3/2:h3:Lx3-h3/2); %Staggered in x2 cell wall x1 x3
[xdz,ydz,zdz] = ndgrid(h1/2:h1:Lx1-h1/2, h2/2:h2:Lx2-h2/2, 0:h3:Lx3); %Staggered in x3 cell wall x1 x2


%% Fields
if mode == 2
    a = 0.01; b = 10;
    u = @(x,y,z)(1 + exp( -( ((x - 0.5).^2)/a + ((y - 0.5).^2)/a + ((z - 0.5).^2)/a) ) );
    up = @(t,x,y,z)( (1 - 2*t) .* (exp( -( ((x - 0.5).^2)/a + ((y - 0.5).^2)/a + ((z - 0.5).^2)/a) )) );
        
    m = @(t,x,y,z)(atan(b*(t - 1/2))./a) .* (exp( -( ((x - 0.5).^2)/a + ((y - 0.5).^2)/a + ((z - 0.5).^2)/a) ));   % Allow for anisotropy with c
    J = @(t,x,y,z)(up(t,x,y,z).*m(t,x,y,z));

end

if mode == 1
    f = sin(pi*xc).*sin(pi*yc).*sin(pi*zc);
    fpp = -3*pi^2*sin(pi*xc).*sin(pi*yc).*sin(pi*zc);
end
%%%%%%%%%%%%%%%%%% OPERATORS
%%%%%%%
% derivatives on walls (going from centre to nodes) <GRAD>
ddxn = @(m,k) 1/k*spdiags([-ones(m+1,1) ones(m+1,1)],[-1,0],m+1,m); 
% derivative function (calculating from nodes to centre) <DIV>
ddxc = @(m,k) 1/k*spdiags([-ones(m+1,1) ones(m+1,1)],[0,1],m,m+1); 
% Average function (Go from centres to walls (will also need ghost))
av = @(m) 0.5*spdiags([ones(m,1) ones(m,1)],[-1,0],m+1,m); 
%%%%%%% 

I   = speye(n);     % Create Identities of Appropriate size
Dn  = ddxn(n,h);  % Create 1D Operators 
Dc  = ddxc(n,h);  % Create 1D Operators 
Av  = av(n);

if mode == 1
    Av([1,end]) = [1,1];
    Dn([1,end]) = [2/h; -2/h];
end
% u Boundary Conditions
% If we change the function to be dependent on y and z will need B2 and B3
% but right now we can just reuse them since we are skipping u(y) and u(z)
% M average boundary conditions 
% m(anistropy constant, up(direction we are differentiating),x,y,z)
% if mode == 2
%     b = sparse(n+1,1);
%     b(1)  = 1/h * u(0);
%     b(end) = -1/h * u(1);
%     
%     c1 = 1; c2 = 2; c3 = 3; % Anisotropy factors
%     mB1 = sparse(n+1,1);   % Create sticks of boundary conditions for m
%     mB2 = mB1;
%     mB3 = mB1;
%     mB1(1)   = 0.25 * 1/m(c1,-h/2,h/2,h/2);
%     mB1(end) = 0.25 * 1/m(c1,1+h/2,1-h/2,1-h/2);
%     mB2(1)   = 0.25 * 1/m(c2,h/2,-h/2,h/2);
%     mB2(end) = 0.25 * 1/m(c2,1-h/2,1+h/2,1-h/2);
%     mB3(1)   = 0.25 * 1/m(c3,h/2,h/2,-h/2);
%     mB3(end) = 0.25 * 1/m(c3,1-h/2,1-h/2,1+h/2);
% end

%%% Ramp up to 2D then to 3D
%% 3D
Dn1 = kron(I,kron(I,Dn));
Dn2 = kron(I,kron(Dn,I));
Dn3 = kron(Dn,kron(I,I));

Dc1 = kron(I,kron(I,Dc));
Dc2 = kron(I,kron(Dc,I));
Dc3 = kron(Dc,kron(I,I));

A1 = kron(I,kron(I,Av));
A2 = kron(I,kron(Av,I));
A3 = kron(Av,kron(I,I));


% Create Large operators
DIV  = [Dc1,Dc2,Dc3];
GRAD = [Dn1;Dn2;Dn3];



if mode == 1
    m1 = ones(n,n,n);
    m1 = m1(:);
    m2 = m1; m3 = m1; 
end
    
    
    
 if mode == 2
 m1 = m(xc,xc,yc,zc);
 m2 = m(yc,xc,yc,zc);
 m3 = m(zc,xc,yc,zc);
 end
 
Am1 = spdiags(A1*m1(:),0,size(A1,1),size(A1,1));
Am2 = spdiags(A2*m2(:),0,size(A2,1),size(A2,1));
Am3 = spdiags(A3*m3(:),0,size(A3,1),size(A3,1));

%Sinv = blkdiag(diag(A1*m1(:)),diag(A2*m2(:)),diag(A3*m3(:)));
Sinv = blkdiag(Am1,Am2,Am3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Solve
% Solving A = DIV*Sinv*GRAD in 2 steps to take advantage of backslash
% operator in matlab
if mode == 1
    SG = Sinv\GRAD;
    A = DIV*SG;
    u = A\fpp(:);

% Residual (computed solution minus function)
    relres = abs(u - f(:))./f(:); % 0.001785
    res = abs(u - f(:));
%u = reshape(u,size(q));

    r(ii) = norm(res,'inf');
end

if mode == 2
    j1 = J(xdx,xdx,ydx,zdx);
    j2 = J(ydy,xdy,ydy,zdy);
    j3 = J(zdz,xdz,ydz,zdz);
    Jvec = [j1(:) ; j2(:) ; j3(:) ];
    %Bvec = [B1; B2; B3];
    U = u(xc,yc,zc);
    r2 = abs(GRAD*U(:) - Sinv*Jvec);
    r(ii) = norm(r2,'inf');
    l = length(r2)*1/3;
    rx = reshape(r2(1:l),n+1,n,n);
    ry = reshape(r2(l+1:2*l),n,n+1,n);
    rz = reshape(r2(2*l + 1: end),n,n,n+1);
    
end
 
 
end

if mode == 1
    figure(1)
    semilogy(nn,r,'--', nn,1./(nn).^2,':', nn, 1./(nn).^3,':', nn,1./(nn).^4,':')
    legend('residual','O(h^2)','O(h^3)','O(h^4)')
    title('Convergence Residual showing O(h^2)')
    ylabel('log residual')
    xlabel('N')
end

if mode == 2
    res = rz;
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
end









