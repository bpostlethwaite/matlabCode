% Example 8.4 -- Figure 8.4 : inverse power method for a simple problem
clear all
close all

u = 1:32; A = diag(u); I = eye(32);

A1 = A-33*I; A2 = A-35*I;

x = ones(32,1); x = x/norm(x);
initial = x'*A*x;

for i = 1:1000, 
    x = A1\x; x = x/norm(x); 
    lamA1(i) = x'*A*x; 
    if abs(lamA1(i)-32)<1e-5 break, end
end

x = ones(32,1); x = x/norm(x);

for i = 1:100, 
    x = A2\x; x = x/norm(x); 
    lamA2(i) = x'*A*x; 
    if abs(lamA2(i)-32)<1e-5 break, end
end

lamA1 = [initial lamA1];
semilogy(0:length(lamA1)-1,abs(lamA1-32),'r--*')
hold on
lamA2=[initial lamA2];
semilogy(0:length(lamA2)-1,abs(lamA2-32),'-*')

legend('\alpha=33','\alpha=35')
xlabel('iteration')
ylabel('absolute error')