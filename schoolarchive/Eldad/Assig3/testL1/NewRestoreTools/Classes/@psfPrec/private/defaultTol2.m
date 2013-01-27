function tol = defaultTol2(E, b)
%
%     tol = defaultTol2(E, b)
%
% Stripped down version of defaultTol.
% This function uses the GCV function to choose a default truncation
% tolerance for constructing a psfPrec.  
%

%  J. Nagy  01-9-03

%e = reshape(E', prod(size(E)), 1);
e=E(:);

[e, idx] = sort(abs(e));
e = flipud(e);

%beta = reshape(fftn(b)', prod(size(b)), 1);
beta=fftn(b);
beta=beta(:);

beta = abs( flipud( beta(idx) ) ) / sqrt(prod(size(b)));

n = length(e);

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
tol = e(reg_min);

