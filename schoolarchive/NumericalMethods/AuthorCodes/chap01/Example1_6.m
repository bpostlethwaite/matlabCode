% Example 1.6 : Evaluate integral using recursive formula

y(1) = log(11) - log(10);  % this is y_0
for n=1:30
    y(n+1) = 1/n - 10*y(n);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For comparison, use numerical quadrature for the same purpose
for n = 1:31
    z(n) = quad(@(x)fun1_6(x,n-1),0,1,1.e-10);
end

format long g
fprintf ('recursion result   quadrature result   abs(difference)\n')
for n = 1:31
    fprintf (' %e       %e      %e\n',y(n),z(n),abs(y(n)-z(n)))
end