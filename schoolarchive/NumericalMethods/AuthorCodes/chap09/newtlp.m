function [dx,dy,ds] = newtlp(A,b,c,x,y,s,mu,dx,ds)
%
% function [dx,dy,ds] = newtlp(A,b,c,x,y,s,mu,dx,ds)
%
% A Newton step for lp
% Use one of three direct elimination methods
%
  
rc = A'*y + s - c;
rb = A*x - b;
rt = x.*s - mu;
if nargin == 9
  rt = rt + dx.*ds;
end

rhs = -rb + A * ((rt - x.*rc)./s);
Mat = A * diag(x./s) * A';
dy = Mat \ rhs;
ds = -rc - A'*dy;
dx = -(x.*ds + rt)./s;