function tol = GCVforSVD2(s, beta)
%
%     tol = GCVforSVD2(s, beta);
%
% Stripped down version of GCVforSVD.
%
% This function uses the GCV function to choose a default truncation
% tolerance for TSVD regularization.
%
% If the blurring matrix has the SVD:  K = U*S*V', then ...
%
% Input:  s - vector containing the singular values of K
%      beta - vector U'*b
%

%  J. Nagy   01-09-03

[s, idx] = sort(abs(s));
s = flipud(s);

beta = abs( flipud( beta(idx) ) ) / sqrt(length(beta));

n = length(s);

%
%  Here we use the GCV function to pick a tolerance.
%
rho2 = zeros(n-1,1);
rho2(n-1) = beta(n)^2;
for k=n-2:-1:1
  rho2(k) = rho2(k+1) + beta(k+1)^2;
end
G = zeros(n-1,1);
for k=1:n-1
  G(k) = rho2(k)/(n - k)^2;
end
[minG,reg_min] = min(G);
tol = s(reg_min);


