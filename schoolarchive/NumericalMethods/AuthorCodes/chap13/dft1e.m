function [a0,a,b] = dft1e(y)
%
%  function [a0,a,b] = dft1e(y)
%
%  Given m values of y, m even, these values are interpreted as 
%  function values at m equidistant points on [0, 2*pi).
%  construct the DFT

 y = y(:); m = length(y); l = m / 2; 
 pi2 = 2*pi; pi2m = pi2/m;
 x = 0: pi2m: (m-1)*pi2m;
 
 a0 = sum(y)/l;
 for k=1:l
    co = cos(k*x); si = sin(k*x);
    a(k) = (co*y)/l;
    if k < l
      b(k) = (si*y)/l;
    end
 end

 
