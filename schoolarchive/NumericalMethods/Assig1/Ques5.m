% Question 5

clear all
close all
clc

x = linspace(-120,120,1000);

fx = x./(cosh(x/4));

plot(x,fx)
title(sprintf('x/cosh(x/4)  note the points where equation \n has one value for given alpha (at peaks and at zero)')) 

f = 'x/4 .* tanh(x/4) - 1';
func = inline(f);


nprobe = 50;
b = 10;
a = -10;
tol = 1e-7;
maxiter = 50;

[roots,iter] = rootgrabber(func,nprobe,a,b,tol,maxiter,true);

x = roots;
fx = x./(cosh(x/4));

for ii = 1:length(fx)
    fprintf('\nfunction has single root: %2.5g at alpha: %2.5g',x(ii),fx(ii))
end
fprintf('\n')

