% Figures 7.5 and 7.7 : comparison of methods for Poisson model problems

clear all
close all
N = 31;
tol = 1.e-6;

% build linear system
A = delsq(numgrid('S',N+2));
n = size(A,1);
b = A*ones(n,1);

% Example 7.10, Figure 7.5

Dv = diag(A);         % diagonal part of A as vector
D = diag(Dv);         % diagonal part of A as matrix
E = tril(A);          % lower triangular part of A
omega = 1+cos(pi/(N+1))/(1+sin(pi/(N+1)))^2; % for SOR
Eomega = ((1-omega)*D + omega*E)/omega;      % for SOR

% Jacobi
x = zeros(n,1); r = b;
for i=1:10000,
  x = x + r./Dv;
  r = b - A*x;
  rJ(i)=norm(r)/norm(b);
  if rJ(i) < tol, break, end
end
xJ = x;

% Gauss Seidel
x = zeros(n,1); r = b;
for i=1:10000,
  x = x + E \ r;
  r = b - A*x;
  rGS(i) = norm(r)/norm(b);
  if rGS(i) < tol,  break, end
end
xGS = x;

% SOR
x = zeros(n,1); r = b;
for i=1:10000,    
  x = x + Eomega \ r;
  r = b - A*x;
  rSOR(i) = norm(r)/norm(b);
  if rSOR(i) < tol, break, end
end
xSOR = x;

% Example 7.11, Figure 7.7

% Conjugate gradient, using Matlab's builtin routines
cpu0 = cputime;
% CG
[xCG,flagCG,relresCG,iterCG,rCG] = pcg(A,b,tol,2000);
cgtime = cputime - cpu0 

% preconditioned CG with incomplete Cholesky IC(0)
R=cholinc(A,'0');
[xCGpc,flagCGpc,relresCGpc,iterCGpc,rCGpc] = pcg(A,b,tol,2000,R',R);
cgic0time = cputime - (cpu0+cgtime)

% preconditioned CG with incomplete Cholesky IC(0.01)
R=cholinc(A,.01);
[xCGpct,flagCGpct,relresCGpct,iterCGpct,rCGpct] = pcg(A,b,tol,2000,R',R);
cgictoltime = cputime - (cpu0+cgtime+cgic0time)

% Plot relative residuals
figure(1)
semilogy(rJ,'m')
hold on
semilogy(rGS,'r-.')
semilogy(rSOR,'--')
legend('Jacobi','Gauss-Seidel','SOR')
xlabel('Iterations')
ylabel('Residual norm')

figure(2)
semilogy(rCG/norm(b),'m')
hold on
semilogy(rCGpc/norm(b),'r-.')
semilogy(rCGpct/norm(b),'b--')
legend('CG','PCG/IC(0)','PCG/IC(.01)')
xlabel('Iterations')
ylabel('Residual norm')