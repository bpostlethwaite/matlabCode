function [x, k, Rnrm, Xnrm, Enrm] = PMR2(A, P, b, x, MaxIts, tol, x_exact)
%
%      [x, k, Rnrm, Xnrm, Enrm] = PMR2(A, P, b, x, MaxIts, tol, x_exact);
%
%  Modified "Preconditioned Minimum Residual Method" for symmetric,
%  possibly indefinte linear systems.
%  Reference: M. Hanke, "Conjugate Gradient Methods for Ill-Posed Problems",
%             Pitman Research Notes in Mathematics, 1995.
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

%  J. Nagy, 04-16-02

%
% check for inputs, and set default values (default tol will be
% set later)
%

disp(' '), disp('Beginning PMR2 iterations')
disp(' '), disp('Iteration number ...')
n = length(b(:));
if nargin < 5
  MaxIts = n;, tol = [];, x_exact = [];
elseif nargin < 6
  tol = [];, x_exact = [];
elseif nargin < 7
  x_exact = [];
end
if isempty(x), x = zeros(size(b));, end
if isempty(MaxIts), MaxIts = n;, end
norm_b = norm(b(:));
if isempty(tol), tol = sqrt(eps)*norm_b;, end

%
%  Initialize some values before iterations begin.
%
Rnrm = zeros(MaxIts, 1);
Xnrm = zeros(MaxIts, 1);

Ax = A*x;

cr = b - Ax;
r = P \ cr;

x_old = x;

v_old = zeros(size(b));
v = P \ (A*r);

u_old = zeros(size(b));
u = A*v;

w_old = zeros(size(b));
w = P \ u;

beta = sqrt(w(:)' * u(:));
v = v / beta;
w = w / beta;
u = u / beta;

if ~isempty(x_exact)
  Enrm = zeros(MaxIts, 1);
  nrm_x_exact = norm(x_exact(:));
  x_error = x - x_exact;
  Enrm(1) = norm(x_error(:)) / nrm_x_exact;
end

Rnrm(1) = norm(cr(:));
Xnrm(1) = norm(x(:));

k = 0;

while ( Rnrm(k+1) > tol & k <= MaxIts-1 )
  k = k + 1;
  fprintf('%5.0f',k)

  rho = r(:)' * u(:);
  x = x + rho*v;
  cr = cr - rho*u;
  r = r - rho*w;

  Rnrm(k+1) = norm(cr(:));
  Xnrm(k+1) = norm(x(:));

  if ~isempty(x_exact)
    x_error = x - x_exact;
    Enrm(k+1) = norm(x_error(:)) / nrm_x_exact;
  end

  Aw = A*w;
  alpha = w(:)'*Aw(:);

  v_new =  w - alpha*v - beta*v_old;
  u_new = Aw - alpha*u - beta*u_old;
  w_new = P\Aw - alpha*w - beta*w_old;

  beta = sqrt(w_new(:)' * u_new(:));
  v_old = v;
  v = v_new/beta;
  w_old = w;
  w = w_new / beta;
  u_old = u;
  u = u_new / beta;

end
Rnrm = Rnrm(1:k+1)/norm_b;
Xnrm = Xnrm(1:k+1);
if ~isempty(x_exact)
  Enrm = Enrm(1:k+1);
else
  Enrm = [];
end