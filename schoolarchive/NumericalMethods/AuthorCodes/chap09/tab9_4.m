% Table 9.4  Using LP program

m = 260; n = 570;
A = randn(m,n);
xh = rand(n,1); c = rand(n,1); b = A*xh;

[x,gap,nbas] = lpm(A,b,c);