% Question 3

clear all
close all

x = 0:.001:1;
f = @(xi)(exp(3*xi).*sin(200*xi.^2) ./ (1 + 20 * xi.^2));
nn = 2.^[4:14];

for jj = 1:length(nn)
    
    n = nn(jj);
    
    % SPLINE COMPUTATION
    xi = 0:1/nn(jj):1;
    PP = spline(xi,f(xi));
    YY = ppval(PP,x);
    err1(jj) = max(abs(f(x) - YY));
      
end



g = 6;
for kk = 1:g
    n = nn(kk);
    a = 0; b = 1;
    ni = 0:1:n;
    ti = cos(pi/(2*(n+1)).* (2*ni(:)+1));
    xc = a + (b-a)/2*(ti+1);
    yc = f(xc);
    
    % CHEBYCHEV LAGRANGE COMPUTATION
    % Lagrange Polynomial Solver
    % Construction
    rho = ones(n+1,1);
    for j=1:n+1
        for i=1:n+1
            if i ~= j
                rho(j) = rho(j) * (xc(j)-xc(i));
            end
        end
    end
    
    w = 1./rho;

    % Evaluation
    for ii = 1:length(x)
        psi = prod( x(ii) - xc);
        p1(ii) =  psi*sum( (w.*yc)./(x(ii) - xc));
    end
    
    err2(kk) = max(abs(p1 - f(x)));
    
end

loglog(nn,err1,':s',nn(1:g),err2,'--^',nn,1./nn.^4)
legend('Cubic Spline max(abs(error))','Chebyshev Lagrange max(abs(error))')
title('loglog of Cubic Spline vs Chebyshev Lagrange interpolation error')
xlabel('Log of Interpolation points')
ylabel('Log max absolute error')