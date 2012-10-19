function c = hermite (t0,t1,f0,f0p,f1,f1p)
%
% use monomial basis to find Hermite cubic
%
A = ones(4,4);
A(1,2) = t0; A(1,3) = t0^2; A(1,4) = t0^3;
A(2,1) = 0; A(2,3) = 2*t0; A(2,4) = 3*t0^2;
A(3,2) = t1; A(3,3) = t1^2; A(3,4) = t1^3;
A(4,1) = 0; A(4,3) = 2*t1; A(4,4) = 3*t1^2;
y = [f0,f0p,f1,f1p]';
c = A \ y;