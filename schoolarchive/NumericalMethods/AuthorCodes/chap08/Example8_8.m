% Example 8.8 : minimum norm solution

format long g
A = [1,0,1; 2,3,5; 5,3,-2;3,5,4; -1,6,3];
b = [4,-2,5,-2,1]';

x = A \ b
res = norm(A*x - b)
normx = norm(x)

AA = [A,sum(A,2)];
xx = AA \ b
resres = norm(AA*xx - b)
normxx = norm(xx)

[U,S,V] = svd(AA);
z = U'*b;
sig = diag(S);
y = z(1:3) ./ sig(1:3);
y = [y;0];
xs = V*y
ress = norm(AA*xs - b)
normxs = norm(xs)

