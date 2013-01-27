function [Q,mesh,fevals] = quads(a, b, tol)
%
% function [Q,mesh,fevals] = quads(a, b, tol)
%
% adaptively evaluate Q - an approximation to the
% the integral from a to b of func(x), within tolerance tol.
%
% mesh is the resulting mesh where f was evaluated
% fevals is the number of function evaluations required.

% initialize

  maxlevel = 50;
  fa = func(a); fb = func(b); fab2 = func((a+b)/2);
  sab = srule (a,b,fa,fab2, fb);

  fevals = 3;
  mesh(1) = a; mesh(2) = b; mesh(3) = (a+b)/2;

% Evaluate the integral

  [Q,mesh,fevals] = quade (a, b, tol, fa, fab2, fb, sab, ...
                           maxlevel, fevals, mesh);
  
% sort mesh in ascending order

  mesh = sort(mesh);
