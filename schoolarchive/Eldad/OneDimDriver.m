%% Variant of Eldad Code for notes in Optmization/Simulation Course
%
% d/dy*J = q;  1/m*J - du/dy = 0
% u(0) = 1;  J(1) = 0    where J = du/dy == flux 
% note u is on staggered mesh, thus no REAL u(0) exists and we must do
% 1/2*(u(ghost) + u((3/2)) = u(0) --> u(ghost) = 2*(u(0) - u(3/2))
% cell size h(k+1/2) = x(k+1) - x(k)

clear all
close all

%% Set constants
kk = 1;

for n = [8,16,32,64,128,256,512]

L = 1; % Our Length Scale (Goes from 0 -> 1)
%n = 8; % Number of nodes
h = L/n; % Size of each cell

%% Set up cells, mesh, Operators
% MESH
yN = 0:h:L; yN=yN(:); % y should be n+1 column vector in length
assert(length(yN) == (n+1),'y is not n+1 in length')
yC = yN(1:end-1) + diff(yN)/2; yC= yC(:);

% OPERATORS
e = ones(n+1,1);
Dn2c = 1/h*spdiags([-e,e],[0,1],n,n+1);
Dc2n = 1/h*spdiags([-e,e],[-1,0],n+1,n);
Ac2n  = 1/2*spdiags([e,e],[-1,0],n+1,n);

%% Fields
a = 10^9; b = pi/2+0.01;
u = @(t)(atan(a*t - 1/2*a).*t - 1/2*atan(a*t - 1/2*a) ...
    -1/(2*a)*log((a*t - 1/2*a).^2 + 1) + b*t);
up = @(t)(atan(a*(t-0.5)) + b);
m = @(t)(1./up(t).*exp(t));
J = @(t)(up(t).*m(t));
rhs = @(t)(exp(t));
%% Boundary Conditions and Boundary Corrections
yB = [yN(1); yN(end)];

s = Ac2n*(1./m(yC));
s(1) = (1/m(h/2) + 1/m(-h/2))/4;
s(end) = (1/m(1-h/2) + 1/m(1+h/2))/4;
Sdiag = diag(s);
SdiagINV = diag(1./s);
% Boundary Condition Matrix
B = sparse(n+1,2); B(1,1) = 1/h; B(end,end) = -1/h;

% check truncation error
r1 = -Dc2n*u(yC) + Sdiag*J(yN) + (B*u(yB));
r2 = Dn2c * J(yN) - rhs(yC);
% now solve and compare solution error
A = Dn2c*SdiagINV*Dc2n;
b = rhs(yC) + Dn2c*SdiagINV*B*u(yB);
uN = A\b; uA = u(yC);
r3 = uA-uN;

%% print some info
rho(kk,:) = [h,norm(r1,'inf'),norm(r2,'inf'),norm(r3,'inf')]; kk = kk+1;
fprintf('%3.2e %3.2e %3.2e %3.2e\n',h,norm(r1,'inf'),norm(r2,'inf'),norm(r3,'inf'));
end
%% check convergence rate
fprintf('\n\n\n Convergence summary\n\n');
fprintf('h |r1|_inf |r2|_inf |error|_inf \n\n')
fprintf('%3.2e %3.2e %3.2e %3.2e\n',(rho(1:end-1,:)./rho(2:end,:))')

figure(1)
t = 1:kk-1;
semilogy(t,rho(:,1),t,rho(:,2),t,rho(:,3),t,rho(:,4))
legend('h','r1','r2','r3','Location','Best')

figure(2)
plot(uN)
