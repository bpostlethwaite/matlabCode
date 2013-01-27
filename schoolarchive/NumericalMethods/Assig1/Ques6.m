clear all
close all
clc

tol = 1e-10;
nmax = 100;
ii = 1;

for n = [8,16,32,64]
    
x = linspace(1,2,n);
fx = log(x)';
y = log(2)*(x-1);

y(1) = 0;
y(end) = log(2);
y = y';

for k = 1:nmax
    [fy,Jy] = buildJaco(y);
    p = -Jy \ fy;
    y(2:end-1) = y(2:end-1) + p;
    plot(x,y)
    if norm(p) < tol*(1+norm(y))
        fy = buildJaco(y);
        break
    end
    
end

E = fx - y;
relE = (fx - y)./fx;

figure(ii)
        subplot(3,1,1)
        plot(x,fx,'--r',x,y)
        legend('ln(x)', 'newton approx','Location','Best')
        title(sprintf(['Real and Newton calculated function for n = %d \n'...
            'and number of iterations: %d'],n,k))
        xlabel('x')
        ylabel('fx')
        
        subplot(3,1,2)
        plot(x,E)
        title(sprintf('Absolute Error for n = %d',n))
        xlabel('x')
        ylabel('Absolute Error')
        
        subplot(3,1,3)
        plot(x,relE)
        title(sprintf('Relative Error for n = %d',n))
        xlabel('x')
        ylabel('relative error')
        
ii = ii + 1;     
    
end