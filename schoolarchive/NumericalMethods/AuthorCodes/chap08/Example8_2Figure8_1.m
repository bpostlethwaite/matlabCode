%Example 8.2 -- Figure 8.1 : Power method for simple examples
clear all
close all

n = 32;
u = 1:n; v = [1:n-2,n-2,n];
A = diag(u); B=diag(v);

lam_exact = n;
x_initial = ones(n,1);

% Do it for A
x = x_initial; x = x/norm(x);
initial = x'*A*x;
for i = 1:1000, 
    x = A*x; 
    x = x/norm(x); 
    lamA(i) = x'*A*x; 
    if abs(lamA(i)-lam_exact)<1e-5 break, end
end

% Do it for B
x = x_initial;
x = x/norm(x);
for i = 1:1000, 
    x = B*x; 
    x = x/norm(x); 
    lamB(i) = x'*B*x; 
    if abs(lamB(i)-lam_exact)<1e-5 break, end
end

lamA = [initial,lamA];
semilogy(0:length(lamA)-1,abs(lamA-lam_exact))
hold on
lamB = [initial,lamB];
semilogy(0:length(lamB)-1,abs(lamB-lam_exact),'r--')
legend('A=diag(u)','B=diag(v)')
xlabel('iteration')
ylabel('absolute error')