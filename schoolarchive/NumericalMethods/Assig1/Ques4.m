clear all
close all
clc

f = 'sinc(x)';
func = inline(f);
nprobe = 50;
b = 10;
a = -10;
tol = 1e-7;
maxiter = 50;

[roots,iter] = rootgrabber(func,nprobe,a,b,tol,maxiter,true);