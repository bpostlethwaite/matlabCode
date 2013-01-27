function x = ainvb_scaled(A,b)
%
% function x = ainvb_scaled(A,b)
%
% solve Ax = b

[p,LU] = plu_scaled (A);
y = forsub (LU,b,p);
x = backsub (LU,y);
