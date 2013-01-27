% Example 5.13 : Cholesky for random matrices

n = 500;
C = randn(n,n); A = C'*C;
xe = randn(n,1); % the exact solution
b = A*xe;        % generate right hand side data

R = chol(A);     % Cholesky factor
% the following line is for compatibility with forsub
D = diag(diag(R)); L = D \ R'; bb = D \ b; p = 1:n;
y = forsub(L,bb,p);            % forward substitution R'y = b
x = backsub(R,y);              % backward substitution Rx = y
rerx = norm(x-xe)/norm(xe)     % error by Cholesky 

xd = ainvb(A,b);               % ignore spd and use partial pivoting
rerxd = norm(xd-xe)/norm(xe)   % error by general routine