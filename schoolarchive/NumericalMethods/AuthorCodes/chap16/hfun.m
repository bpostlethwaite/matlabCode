function H = hfun(t,y)
%
%  Hamiltonian invariant H = hfun(t,y)

isw = 2;

if isw == 1
   %  Henon-Heiles

   T = .5*(y(3)^2+y(4)^2);
   V1 = .5*(y(1)^2+y(2)^2);
   V2 = y(1)^2*y(2) - y(2)^3/3;
   H = T + V1 + V2;

elseif isw == 2
    
    %  FPU
    m = 3; omega = 100;
    T = .5* norm(y(2*m+1:4*m))^2;
    V1 = .5* omega^2 * norm(y(m+1:2*m))^2;
    V2 = .25*( (y(1)-y(m+1))^4 + (y(m)+y(2*m))^4);
    V3 = .25*sum( (y(2:m)-y(1:m-1)-y(m+2:2*m)-y(m+1:2*m-1)).^4 );
    H = T + V1 + V2 + V3;
end
