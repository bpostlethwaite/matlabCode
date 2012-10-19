% Question 3b
% Simple diffusion discretizatoin
clear all
close all
n = 100;
h = 1/n;

l = -ones(n-1,1)*1/h^2;
l(end) = l(end)*2;
c = 2*ones(n,1)*1/h^2;
r = -ones(n-1,1)*1/h^2;

t = h:h:1;
b = (pi/2)^2 * sin(pi/2 *t);

v = tridiagSolv(l,c,r,b);

u = sin(pi/2*t);

TITLE = sprintf('\nInfinite norm of ||v-u|| is %2.3g\n',max(abs(v-u)));
TITLE

plot(t,v,'p',t,u,'r^')
legend('calculated','actual')
title(TITLE)
