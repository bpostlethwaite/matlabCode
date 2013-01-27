clear all
close all
clc

f = '2*cosh(x/4) - x';
func = inline(f);
nprobe = 50;
b = 10;
a = -10;
tol = 1e-7;
maxiter = 50;

[roots,iter] = rootgrabber(func,nprobe,a,b,tol,maxiter,true);


