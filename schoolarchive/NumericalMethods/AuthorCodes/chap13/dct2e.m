function yy = dct2e(xx,a,method)
%
%  function yy = dct2e(xx,a,method)
%
%  Given DCT coefficients a, evaluate triginometric polynomial at xx.

   m = length(a); n = m-1;
   if method == 1   % our way
      yy = .5 * (a(1)*ones(size(xx)) + a(m)*cos(n*xx)); 
      for k=1:n-1
         yy = yy + a(k+1)*cos(k*xx);
      end
      
   elseif method == 2  % DCT-II (III)
     yy = .5 * a(1)*ones(size(xx)); 
      for k=1:n
         yy = yy + a(k+1)*cos(k*xx);
      end
   end



 