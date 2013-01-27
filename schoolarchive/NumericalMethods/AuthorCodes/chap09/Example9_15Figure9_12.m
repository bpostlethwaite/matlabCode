% Example 9.15 -- Figure 9.12 : underdetermined solutions minimizing 
% l_2 and l_1 solution norms. The one using l_1 is sparse.

clear all
close all
% create data
m = 260; n = 570;
J = randn(m,n);

xh = zeros(n,1);
for i = 20:20:260
    xh(i) = 1 - 10/i;
end
b = J*xh;

% minimum l_2 norm solution
[U,S,V] = svd(J);
bb = U'*b; z = zeros(n,1);
for i = 1:m
  z(i) = bb(i) / S(i,i);
end
y2 = V*z;
res = norm(b - J*y2)
y2n2 = norm(y2)
y2n1 = norm(y2,1)
plot (1:n,y2')
hold on

% minimum l_1 norm solution
c = ones(2*n,1);
A = .5*[J,-J];
kk = 20;

[w,gap,nbas] = lpl1 (A,b,c,kk);

y1 = .5*(w(1:n) - w(n+1:2*n));
res = norm(b - J*y1)
y1n2 = norm(y1)
y1n1 = norm(y1,1)
plot (1:n,y1','g')
legend('l_2 minimizer','l_1 minimizer')
xlabel('component i')
ylabel('solution y')
