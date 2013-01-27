function tol = GCVTolPlot(A, b)
%
%        tol = GCVTolPlot(A, b)
%
% Given a psfMatrix A, and RHS image b, this function uses the
% generalized cross validation (GCV) method to help find a tolerance
% for the preconditioner.
%
% The GCV function is:
%          || P*x - b ||^2
%    G = -------------------
%        (trace(I - P*P_I)^2
% where, P_I is the truncated inverse of the circulant
% preconditioner, P.

% Jim Nagy, 12/30/01
% This code is a modification of P.C. Hansen's gcv function
% in regularization tools.

% Initialization.
P = psfPrec(A, b, 0);
e = 1 ./ P.matdata;
e = reshape(e', prod(size(e)), 1);
[e, idx] = sort(abs(e));
e = flipud(e);
beta = abs( reshape(fft2(b)', prod(size(b)), 1) );

n = length(e);

rho2 = zeros(n-1,1);

rho2(n-1) = beta(n)^2;

for k=n-2:-1:1
  rho2(k) = rho2(k+1) + beta(k+1)^2; 
end

G = zeros(n-1,1);
for k=1:n-1
  G(k) = rho2(k)/(n - k)^2;
end

figure
reg_param = [1:n-1]';
semilogy(reg_param,G), xlabel('k'), ylabel('G(k)')

[minG,reg_min] = min(G);
tol = e(reg_min);

hold on
semilogy(reg_min,minG,'o',[reg_min,reg_min],[minG/1000,minG],'--')
legend('GCV function', sprintf('tol = %f', tol))
hold
title('GCV function')




