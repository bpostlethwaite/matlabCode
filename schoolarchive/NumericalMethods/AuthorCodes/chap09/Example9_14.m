% Example 9.14 : a larger lp problem. Note that because random numbers are
% used to generate data this will not reproduce Table 9.3.

m = 260; n = 570;
A = randn(m,n);
xh = rand(n,1); c = rand(n,1); b = A*xh;

[x,gap,nbas] = lpm (A,b,c);
