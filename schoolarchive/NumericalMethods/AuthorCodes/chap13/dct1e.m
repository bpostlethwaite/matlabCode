function [a] = dct1e(y,method)
%
%  function [a] = dct1e(y,method)
%
%  Given m values of y, m even, these values are interpreted as 
%  function values at m equidistant points on [0, 2*pi).
%  construct the DCT

 y = y(:); m = length(y); n = m-1;
 
 if method == 1  % our way
   % m odd, both ends included  
   pi2 = 2*pi; pi2m = pi2/n;
   x = 0: pi2m: n*pi2m;
 
   for k=0:n
     co = cos(k*x);
     a(k+1) = co*y;
   end
   a = a / n;
   
 elseif method == 2  %  dct-II from wikipedia
   pim = pi/m;
   x = .5*pim: pim: (n+.5)*pim;
 
   for k=0:n
     co = cos(k*x);
     a(k+1) = co*y;
   end
   a = a / m*2;
   
 end
 
