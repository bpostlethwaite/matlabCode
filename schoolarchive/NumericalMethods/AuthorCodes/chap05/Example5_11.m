% Example 5.11 : comparing scaled and unscaled GEPP for a special example

n = 100; h = 1/(n-1); K = 100;
A = zeros(n,n);
for i = 2:n-1
  A(i,i) = -2/h^2 - K;
  A(i,i-1) = -1/h^2; A(i,i+1) = -1/h^2;
end
A(1,1) = 1; A(n,n) = 1;  % end definition of A

xe = ones(n,1);          % exact solution of 1's  
b = A*xe;                % corresponding right hand side

%solve using ainvb
xu = ainvb(A,b);
err_ainvb = norm(xu-xe)

%solve using scaled GEPP version ainvb_scaled
xs = ainvb_scaled(A,b);
err_ainvb_scaled = norm(xs-xe)