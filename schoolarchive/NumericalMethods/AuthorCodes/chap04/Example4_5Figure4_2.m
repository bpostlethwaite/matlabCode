% Figure 4.2

A = zeros(5,5);
A(1,2) = 2; A(2,1) = -2; A(3,3) = -1; A(3,4) = 2.5; 
A(4,3) = -2.5; A(4,4) = -1; A(5,5) = -2;

B = [1 2 3;4 5 6;7 8 9;-1 -2 -3;-4 -5 -6];
[Q,R] = qr(B);

A = Q'*A*Q
lam = eig(A);
plot(real(lam),imag(lam),'ro')
hold on
xlabel('\Re(\lambda)')
ylabel('\Im(\lambda)')
axis([-5 1 -5 5])
plot([-5 1],[0 0],'k')
plot([0 0], [-5 5], 'k')