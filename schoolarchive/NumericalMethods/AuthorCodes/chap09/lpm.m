function [x,gap,nbas] = lpm (A,b,c)
%
% function [x,gap,nbas] = lpm (A,b,c)
%
% solve the linear programming problem
% min c^T x  s.t. Ax = b, x >= 0 . 
%
% A is l x m,  b is l x 1, c is m x 1.
% return solution x and duality gap 
% (should be close to 0 if all well).
% Also, nbas is the number of hops to check basic solutions
% before optimality is reached.

[l,m] = size(A);
scaleb = norm(b) + norm(A,inf) + 1;
scalec = norm(c) + norm(A,inf) + 1;
tolf = 0.01; otol = 1-tolf; 
toln = 1.e-9; tolp = 1.e-10; tolb = 1.e-4;
nbas = 0;

% Initial guess
x = ones(m,1); s = ones(m,1); y = zeros(l,1);

fprintf('itn       gap          infeas            mu\n')
for it = 1:2*m+10 % iterate, counting pred-cor as one

  % duality measure
  mu = (x'*s)/m;
  % predict correction
  [dx,dy,ds] = newtlp(A,b,c,x,y,s,0);

  % incorporate positivity into constraints
  alfa = 1; beta = 1;
  for i=1:m
    if dx(i) < 0, alfa = min(alfa, -x(i)/dx(i)); end
    if ds(i) < 0, beta = min(beta, -s(i)/ds(i)); end
  end

  % the would-be duality measure
  muaff = ( (x+alfa*dx)' * (s+beta*ds) ) / m;
  % centering parameter
  sigma = (muaff/mu)^3;

  % correct towards center path
  smu = sigma * mu;
  [dx,dy,ds] = newtlp(A,b,c,x,y,s,smu,dx,ds);

  % incorporate positivity into constraints
  alfa = 1; beta = 1;
  for i=1:m
    if dx(i) < 0, alfa = min(alfa, -otol*x(i)/dx(i)); end
    if ds(i) < 0, beta = min(beta, -otol*s(i)/ds(i)); end
  end

  % update solution
  x = x + alfa*dx;
  s = s + beta*ds;
  y = y + beta*dy;

  % check progress
  infeas = norm(b - A*x)/scaleb + norm(c - A'*y - s)/scalec;
  gap = (c'*x - b'*y) / m;
  if (infeas > 1.e+12)+(gap < -toln)
    fprintf('no convergence: perhaps no solution')
    return
  end 
  fprintf('%d    %e    %e    %e\n',it,gap,infeas,mu)

  if (abs(infeas) < toln)*(abs(gap) < toln), return, end

  % hop to next basic solution
  if gap < tolb
    nbas = nbas + 1;
    [xx,sortof] = sort(-x);
    [xx,yy,ss] = basln(A,b,c,sortof);
    gap = (c'*xx - b'*yy) / m;
    if (sum(xx+tolp >= 0) > m-1)*(sum(ss+tolp >= 0) > m-1)...
            *(abs(gap) < toln)
      x = xx; 
      return
    end 
  end

end