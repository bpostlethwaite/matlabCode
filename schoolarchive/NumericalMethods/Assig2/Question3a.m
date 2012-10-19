%% QUestion 3
clear all
close all

n = 100000;

l = -1:-1:-n;
r = -1:-1:-n;
c = 3:3:n*3;
b = ones(n,1);

x = tridiagSolv(l(1:n-1),c,r(2:n),b);
