%3D Poisson Driver
clear all
close all

%% Set constants
nn = 10;
L = 1;

n1 = nn; % Number of cells in x1 direction
n2 = nn; % Number of cells in x2 direction
n3 = nn; % Number of cells in x3 direction

h1 = L/n1; % cell length in x1 direction
h2 = L/n2;
h3 = L/n3;

Lx1 = n1*h1;
Lx2 = n2*h2;
Lx3 = n3*h3;

% Useful 3D shortcuts
h = h1*h2*h3;
n = n1*n2*n3;
e  = ones(n,1);

%% Grid
% MESH
[x,y,z] = ndgrid(0:h1:Lx1,0:h2:Lx2,0:h3:Lx3); % Cell nodes
[xc,yc,zc] = ndgrid(h1/2:h1:Lx1-h1/2, h2/2:h2:Lx2-h2/2, h3/2:h3:Lx3-h3/2); % Cell centres
[xdx,ydx,zdx] = ndgrid(0:h1:Lx1, h2/2:h2:Lx2-h2/2, h3/2:h3:Lx3-h3/2); %Staggered in x1 cell wall x2 x3
[xdy,ydy,zdy] = ndgrid(h1/2:h1:Lx1-h1/2, 0:h2:Lx2, h3/2:h3:Lx3-h3/2); %Staggered in x2 cell wall x1 x3
[xdz,ydz,zdz] = ndgrid(h1/2:h1:Lx1-h1/2, h2/2:h2:Lx2-h2/2, 0:h3:Lx3); %Staggered in x3 cell wall x1 x2

%% 1D OPERATORS
% derivatives on walls (going from centre to nodes) <GRAD>
ddxn = @(m,k) 1/k*spdiags([-ones(m+1,1) ones(m+1,1)],[-1,0],m+1,m);
% derivative function (calculating from nodes to centre) <DIV>
ddxc = @(m,k) 1/k*spdiags([-ones(m+1,1) ones(m+1,1)],[0,1],m,m+1);
% Average function (Go from centres to walls (will also need ghost))
av = @(m) 0.5*spdiags([ones(m,1) ones(m,1)],[-1,0],m+1,m);
av2 =  @(m) spdiags([ones(m,1) ones(m,1)],[0,1],m,m+1);

I1   = speye(n1);     % Create Identities of Appropriate size
I2   = speye(n2);     % Create Identities of Appropriate size
I3   = speye(n3);     % Create Identities of Appropriate size

Dn1  = ddxn(n1,h1);  % Create 1D Operators
Dn2  = ddxn(n2,h2);  % Create 1D Operators
Dn3  = ddxn(n3,h3);  % Create 1D Operators

Gn1  = ddxn(n1,h1);  % Create 1D Operators
Gn2  = ddxn(n2,h2);  % Create 1D Operators
Gn3  = ddxn(n3,h3);  % Create 1D Operators

% Boundary Conditions on Grad
Dn1([1,end]) = [2/h1; -2/h1];
Dn2([1,end]) = [2/h2; -2/h2];
Dn3([1,end]) = [2/h3; -2/h3];

Gn1([1,end]) = [0; 0];
Gn2([1,end]) = [0; 0];
Gn3([1,end]) = [0; 0];

Dc1  = ddxc(n1,h1);
Dc2  = ddxc(n2,h2);
Dc3  = ddxc(n3,h3);

Av1  = av(n1);
Av2  = av(n2);
Av3  = av(n3);

Av21  = av2(n1);
Av22  = av2(n2);
Av23  = av2(n3);

Av1([1,end]) = [1,1];
Av2([1,end]) = [1,1];
Av3([1,end]) = [1,1];

%% 3D Operators
DN1 = kron(I3,kron(I2,Dn1));
DN2 = kron(I3,kron(Dn2,I1));
DN3 = kron(Dn3,kron(I2,I1));

GN1 = kron(I3,kron(I2,Gn1));
GN2 = kron(I3,kron(Gn2,I1));
GN3 = kron(Gn3,kron(I2,I1));

DC1 = kron(I3,kron(I2,Dc1));
DC2 = kron(I3,kron(Dc2,I1));
DC3 = kron(Dc3,kron(I2,I1));

A1 = kron(I3,kron(I2,Av1));
A2 = kron(I3,kron(Av2,I1));
A3 = kron(Av3,kron(I2,I1));

A21 = kron(I3,kron(I2,Av21));
A22 = kron(I3,kron(Av22,I1));
A23 = kron(Av23,kron(I2,I1));

AV = [A1;A2;A3];
AV2 = [A21,A22,A23];


D  = [DC1,DC2,DC3];
G  = [DN1;DN2;DN3];
G2 = [GN1;GN2;GN3];

%%  Define q and Q


% Define numerous transmitter electrode locations q
line = 3; % y direction on slices
dipole = [4,6]; % x direction on slices
spacing = [1,2,3,5];
numq = length(line)*length(dipole)*length(spacing);
qtot = zeros(n,numq);
kkk = 1;
for i = 1:length(line);
    for j = 1:length(dipole)
        for k = 1:length(spacing)
            q = zeros(n1,n2,n3);
            q(line(i),dipole(j),1) = -1/(h1*h2*h3);
            q(line(i) + spacing(k),dipole(j),1) = 1/(h1*h2*h3);
            qtot(:,kkk) = q(:);
            kkk=kkk+1;
        end
    end
end
q = qtot(:);

% Define matrix Q which outlines where data will be collected
% If data is collected at every surface cell:
% Cut is cutoff, how many cells from the edge should the data be cutoff at
cut = 1;
cut = cut*2;
Qx = spdiags(ones(n1-cut,1),cut/2,n1-cut,n1);
Qy = spdiags(ones(n2-cut,1),cut/2,n2-cut,n2);
Q = kron(Qy,Qx);
%Q = speye(n1*n2); % Cut off zero.
k = zeros(1,n3); k(1) = 1;
Q = kron(k,Q);
Q = kron(speye(numq),Q); % For multiple q purposes

%% Define model m, add noise
%load model.mat

%
noiselevel = 0.2;
m = ones(size(xc));
mref = m(:);
hr = ceil(nn/2);
m(hr - 1 : hr + 1 ,hr - 1 : hr + 1 , 3:4 ) = 5;
m  = m(:);
mnoise = m(:) + noiselevel*randn(n,1);
%}
%% Matrix Equations

Dmq = kron(speye(numq),D);
Gmq = kron(speye(numq),G);
Avmq = kron(ones(numq,1),AV);
dl = size(Avmq,1);

A       = @(m,A1,A2,A3)(D*DiagAvgM(m,A1,A2,A3)*G);
Amq     = @(m,A1,A2,A3)(kron(speye(numq),A(m,A1,A2,A3)));
S_inmq  = @(m,A1,A2,A3)(kron(speye(numq),DiagAvgM(m,A1,A2,A3)));
dCdu    = @(u,m)(Dmq*S_inmq(m,A1,A2,A3)*Gmq);
diag_m  = @(m)(spdiags(1./m.^2,0,n,n));
dCdm    = @(u,m,A1,A2,A3)(Dmq*(spdiags(Gmq*u,0,dl,dl))*S_inmq(m,A1,A2,A3).^2*Avmq* diag_m(m) );
Um      = @(m,q,A1,A2,A3)( Amq(m,A1,A2,A3) \ q(:) );
Rm      = @(m,e,h)( 1/2 * h * e' * AV2 * ((G2*m).^2));
gradRm  = @(m,e,h)(G2' * spdiags( (AV2' * h * e),0,size(G2,1),size(G2,1)) * G2 *m);
grad2Rm = @(h,e)(G2' * spdiags( (AV2' * h * e),0,size(G2,1),size(G2,1)) * G2);
J       = @(u,m,Q,A1,A2,A3) (-Q * (dCdu(u(:),m)\dCdm(u(:),m,A1,A2,A3)) );

%% Compute datasets and noise


dclean = Um(m,q,A1,A2,A3);
d = Um(mnoise,q,A1,A2,A3);
d = Q*d;
dclean = Q*dclean;
ddirty = 0.05*dclean.*randn(length(dclean(:)),1);
ddirty = dclean + ddirty;
stdevD =  std(ddirty - dclean)
stdev = std(d - dclean)
%dev = noiselevel;% * d(:);
figure(1235)
    plot(ddirty - dclean,'r')
    hold on
    plot(d - dclean)
%noise = dev .* randn(length(d(:)),1);
%d = d(:) + noise;
%d = d(:);

%% View source locations q (change column of qtot to see diff q locations)
% data

data = reshape(d,(n1-cut)*(n2-cut),numq);
qplot = 0;
dplot = 0;


for i = 1:numq;
    qplot = reshape(qtot(:,i),n1,n2,n3);
    figure(546)
    subplot(2,4,i)
        imagesc(qplot(:,:,1))  
        title('Dipole Setup')
        xlabel('X direction')
        ylabel('Y direction')
end




%}

%% Optimize Alpha & Tikhonov curve

% Alpha Find and Tikhonov Curve
%{

alpha = logspace(-8,1,20);

%while  A > Sg + Sg*tol || A < Sg - Sg*tol
for ii = 1:length(alpha)
    mit = mref;
    iter = 1;
    fprintf('ITERATION NUMBER %i\n',ii)
    while  iter <= 4
        
        fprintf('Gauss Newton Step %i\n',iter)
        u = Um(mit,q,A1,A2,A3);
        Y = J(u,mit,Q,A1,A2,A3)' * (Q * Um(mit,q,A1,A2,A3) - d) + alpha(ii) * gradRm(mit - mref,e,h);
        H = J(u,mit,Q,A1,A2,A3)' * J(u,mit,Q,A1,A2,A3) + alpha(ii) * grad2Rm(h,e);
        s = -H\Y;
        mit = mit + s;
        iter = iter + 1;
        fprintf('norm of s is %f\n',norm(s))
        
    end
    
    phid = (Q * Um(mit,q,A1,A2,A3) - d) / stdev;
    phiD(ii) = phid' * phid;
    phiM(ii) = Rm(mit,e,h);
    Mit(:,ii) = mit;
    
    
end

for ii = 1 : length(alpha)
    mnorm(ii) = norm(m - Mit(:,ii));
end

save('Vars.mat','phiD','phiM','Mit','alpha','mnorm')

figure(343)
plot(phiM(2:end),phiD(2:end))
xlabel('Regularization')
ylabel('Data Misfit')
title('Tikhonov Curve')

%}
%% Inversion with optimal Alpha
%
alpha =  7.8e-7;
iter = 1;
mit = mref;
s = 10;
while  iter <= 8
    
    fprintf('Gauss Newton Step %i\n',iter)
    u = Um(mit,q,A1,A2,A3);
    Y = J(u,mit,Q,A1,A2,A3)' * (Q * Um(mit,q,A1,A2,A3) - d) + alpha * gradRm(mit - mref ,e,h);
    H = J(u,mit,Q,A1,A2,A3)' * J(u,mit,Q,A1,A2,A3) + alpha * grad2Rm(h,e);
    s = -H\Y;
    figure(4)
    plot(s)
    mit = mit + s;
    iter = iter + 1;
    snorm(iter) = norm(s);
    mnorm(iter) = norm(m - mit)
    fprintf('norm of s is %f\n',snorm(iter))
end

phid = (Q * Um(mit,q,A1,A2,A3) - d) / stdev;
phiD = phid' * phid
%}


%% Plot
%model = m(:) + 0.001 * randn(n1*n2*n3,1);

model = mnoise;
model = reshape(model,n1,n2,n3);
mInv = reshape(mit,n1,n2,n3);

figure(231);
for ii = 1:9
    subP = subplot(3,3,ii);
    %contourf(m(:,:,ii))
    imagesc(model(:,:,ii))
    title(sprintf('Depth z = %d',(ii-1)))
    caxis([min(min(min(model))) max(max(max(m)))])
    pos=get(subP, 'Position');
    set(subP, 'Position', [pos(1)+0.05 pos(2) pos(3) pos(4)])
    %cax = caxis; %Turn on for absolute color cross over between fig 1
    %and fig 2 (need to turn on corrisponding caxis in fig 2)
end

B=colorbar ('FontSize',12);
set(B, 'Position', [.0314 .11 .0581 .8150])
ylabel(B,'Conductivity anomaly [dimensionless]')

figure(242)
for ii = 1:9
    subP = subplot(3,3,ii);
    %contourf(mInv(:,:,ii))
    imagesc(mInv(:,:,ii))
    title(sprintf('Depth z = %d',(ii-1)*10))
    caxis([min(min(min(mInv))) max(max(max(mInv)))])
    pos=get(subP, 'Position');
    set(subP, 'Position', [pos(1)+0.05 pos(2) pos(3) pos(4)])
    %caxis(cax)  %Turn on for absolute color cross over between fig 1
    %and fig 2
end
B=colorbar ('FontSize',12);
set(B, 'Position', [.0314 .11 .0581 .8150])
ylabel(B,'Conductivity anomaly [dimensionless]')
%}
