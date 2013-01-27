function [x, k, Rnrm, Xnrm, Enrm] = CGLS(A, b, x, MaxIts, tol, x_exact)
%
%      [x, k, Rnrm, Xnrm, Enrm] = CGLS(A, b, x, MaxIts, tol, x_exact);
%
%  Conjugate Gradient for Least Squares
%  Reference: A. Bjorck, "Numerical Methods for Least Squares Problems"
%             SIAM, 1996, pg. 289.
%
%  Input: A  -  object defining the coefficient matrix.
%         b  -  Right hand side vector
%         x  -  initial guess 
%
%  Optional Intputs:
%    MaxIts  -  number of iterations to perform (default = length(b(:)))
%       tol  -  tolerance for stopping (default = sqrt(eps)*norm(b))
%   x_exact  -  if the exact solution is known, we can compute relative
%               errors at each iteration
%
%  Output:
%         x  -  solution
%         k  -  actual number of iterations performed
%      Rnrm  -  norm of the residual at each iteration
%      Xnrm  -  norm of the solution at each iteration
%      Enrm  -  norm of the true error at each iteration (if x_exact is known)
%

%
% J. Nagy, 02-10-01

%
% check for inputs, and set default values (default tol will be
% set later)
%

%disp(' '), disp('Beginning CGLS iterations')
%disp(' '), disp('Iteration number ...') 
h = waitbar(0, 'Beginning CGLS iteration ...');
n = length(b(:));
if nargin < 4
  MaxIts = n;, tol = [];, x_exact = [];
elseif nargin < 5
  tol = [];, x_exact = [];
elseif nargin < 6
  x_exact = [];
end
if isempty(x), x = zeros(size(b));, end
if isempty(MaxIts), MaxIts = n;, end

%
%  Initialize some values before iterations begin.
%
Rnrm = zeros(MaxIts, 1);
Xnrm = zeros(MaxIts, 1);

trAb = A'*b;

nrm_trAb = norm(trAb(:));

if isempty(tol) 
  tol = sqrt(eps)*nrm_trAb;
end

if ~isempty(x_exact)
  Enrm = zeros(MaxIts, 1);
  nrm_x_exact = norm(x_exact(:));
  x_error = x - x_exact;
  Enrm(1) = norm(x_error(:)) / nrm_x_exact;
end

s = A*x;

s = b - s;
r = A'*s;

gamma = r(:)'*r(:);

Rnrm(1) = sqrt(gamma);
Xnrm(1) = norm(x(:));
k = 0;
while ( Rnrm(k+1) > tol & k <= MaxIts-1 )
  k = k + 1;
%  fprintf('%5.0f',k)
  waitbar(k/MaxIts, h)
  if ( k == 1 )
    p = r;
  else 
    beta = gamma / oldgamma;
    p = r + beta * p;
  end
  q = A*p;
  nq = q(:)' * q(:);
  alpha = gamma / nq;
  x = x + alpha * p;
  if ~isempty(x_exact)
    x_error = x - x_exact;
    Enrm(k+1) = norm(x_error(:)) / nrm_x_exact;
  end
  s = s - alpha * q;
  r = A'*s;
  oldgamma = gamma;
  gamma = r(:)' * r(:);
  Rnrm(k+1) = sqrt(gamma);
  Xnrm(k+1) = norm(x(:));
end
Rnrm = Rnrm(1:k+1)/nrm_trAb;
Xnrm = Xnrm(1:k+1);
if ~isempty(x_exact)
  Enrm = Enrm(1:k+1);
else
  Enrm = [];
end
close(h)
