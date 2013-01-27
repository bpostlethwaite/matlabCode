function [xc, yc, ycp, D] = chebdif (f,a,b,n,meth)
%
% function [xc, yc, ycp, D] = chebdif (f,a,b,n,meth)
%
% given a function f on an interval [a,b], find interpolating
% polynomial of degree n at 
% either Chebyshev points (meth=1) or Chebyshev extremum pts (meth=2)
% and differentiate. Use Lagrange form.
% D = differentiation matrix.

% Construct abscissae and data
ni = 0:1:n;
if meth == 1
  ti = cos(pi/(2*(n+1)).* (2*ni+1)); ti = ti(:);
else
  ti = cos(pi/n.* ni); ti = ti(:);
end
xc = a + (b-a)/2*(ti+1);
yc = f(xc);

% Coefficients for Lagrange form
D = zeros(n+1,n+1);
for i=1:n+1
    for j=1:n+1
	if i ~= j, D(i,j) = 1/(xc(j)-xc(i)); end
        for k=1:n+1
          if i == j
	        if k ~= j, D(j,j) = D(j,j) + 1/(xc(j)-xc(k)); end
          else
            if (k ~= j)*(k ~= i)
               D(i,j) = D(i,j)*(xc(i)-xc(k))/(xc(j)-xc(k)); 
            end
          end
        end
    end
    ycp(i) = D(i,:)*yc;
end
ycp = ycp(:);
