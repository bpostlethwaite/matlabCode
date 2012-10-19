% Question 2 Interpolate the Runge Function

clear all
close all


x = -1:.001:1;
f = @(xi)(1./(1 + 25*xi.^2));
nn = 10:10:170;

for jj = 1:length(nn)
    n = nn(jj);
    
    ni = 0:1:n;
    xCheb = ( cos(pi/(2*(n+1)).* (2*ni+1)) )';
    yCheb = f(xCheb) ;
    
    % Lagrange Polynomial Solver
    % Construction
    rho = ones(n+1,1);
    for j=1:n+1
        for i=1:n+1
            if i ~= j
                rho(j) = rho(j) * (xCheb(j)-xCheb(i));
            end
        end
    end
    
    w = 1./rho;

    % Evaluation
    for ii = 1:length(x)
        psi = prod( x(ii) - xCheb);
        p1(ii) =  psi*sum( (w.*yCheb)./(x(ii) - xCheb));
    end
    
    err(jj) = max(abs(p1 - f(x)));
    
    
    
end

% Spectral Accuracy: determine q
y = log(err(:));
A = [ones(length(nn),1),-nn(:)];
c = A\y;
q = exp(c(2));

figure(2)
    plot(x,f(x),':',xCheb,yCheb,'*',x,p1)
    title('Function, Polynomial Interpolation and Chebyshev points')
    xlabel('x axis')
    ylabel('f(x)')
    legend('f(x)','Chebychev Points','Lagrange Interp')
figure(3)
    semilogy(nn,err,nn,exp(c(1)) .* q.^(-nn),'s')
    legend('max error',sprintf('O(q^-^n) with q = %1.2f',q))
    title('semilogy plot of error and spectral accuracy')
    xlabel('number of interpolants')
    ylabel('max absolute error')

