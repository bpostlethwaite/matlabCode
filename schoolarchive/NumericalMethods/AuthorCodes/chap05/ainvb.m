function x = ainvb(A,b)
%
% function x = ainvb(A,b)
%
% solve Ax = b

[p,LU] = plu (A);
y = forsub (LU,b,p);
x = backsub (LU,y);
