% Example 9.7 -- Figure 9.6 : comparing SD, LSD and CG for Poisson's equation.
clear all
close all

N = 31;
tol = 1.e-6;

A = delsq(numgrid('S',N+2));
n = size(A,1);
b = A*ones(n,1);

% Steepest descent
x = zeros(n,1); r = b;
for i=1:10000,
    s = A*r;
    alfa = (r'*r) / (r'*s);
    x = x + alfa*r;
    r = b - A*x;
    rSD(i)=norm(r)/norm(b);
    if rSD(i) < sqrt(tol), break, end
end
xSD = x;

% LSD: gradient descent a la Barzilai-Borwein
x = zeros(n,1); r = b; alfa = (r'*r) / (r'*(A*r));
for i=1:10000,
    alfo = alfa;
    s = A*r;
    alfa = (r'*r) / (r'*s);
    x = x + alfo*r;
    r = b - A*x;
    rLSD(i)=norm(r)/norm(b);
    if rLSD(i) < tol, break, end
end
xLSD = x;

% Conjugate gradient, using Matlab's builtin routines

% CG
[xCG,flagCG,relresCG,iterCG,rCG] = pcg(A,b,tol,2000);

% Plot relative residuals
semilogy(rSD,'m')
hold on
semilogy(rLSD,'r-.')
semilogy(rCG,'--')
legend('SD','LSD','CG')
xlabel('Iterations')
ylabel('Residual')
