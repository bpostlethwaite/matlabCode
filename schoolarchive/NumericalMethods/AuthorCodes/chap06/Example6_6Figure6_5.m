% Example 6.6 : Check execution time of ``backslash'' for ls against 
% ``backslash'' on normal eqns for overdetermined systems

for n = 300:100:1000
  % fill a rectangular matrix A and a vector b with random numbers
  % hoping that A'*A is nonsingular
  m = n+1; % or m = 3*n+1 or something else
  A = randn(m,n); b = randn(m,1);
    
  % solve and find execution times. First, matlab way using QR
  t0 = cputime;
  xqr = A \ b;
  temp = cputime;
  tqr(n/100-2) = temp - t0;
    
  % next use normal equations
  t0 = temp;
  B = A'*A; y = A'*b;
  xne = B \ y;
  temp = cputime;
  tne(n/100-2) = temp - t0; 
end

ratio = tqr./tne;
plot(300:100:1000,ratio)
hold on
xlabel('n')
ylabel('CPU time ratio')