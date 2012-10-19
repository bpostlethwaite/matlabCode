% Example 9.5 -- Table 9.2

% exact solution
xe = [3,.5]';
phie = phi(xe);
nmax = 20; tol = 1.e-8;

% first initial guess
x = [8,.2]';
% second initial guess
x = [8,.8]';

fprintf ('k      ||x_k - x*||      phi_k - phi*        f_kp_k  \n')

% apply Newton's method and print table row after each step
for k=1:nmax
  [fx,Jx] = feval(@funv,x);
  p = -Jx \ fx;
  fprintf ('%d     %e     %e     %e \n',k-1,norm(x-xe),phi(x)-phie, fx'*p )
  x = x + p;  
  if norm(p) < tol*(1+norm(x)) break, end
end
