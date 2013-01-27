clear all; close all;

%{
n = 6;
N = 5000;

x = linspace(-2,n,N);
fx = sin(20*x) + sin(40*x) + sin(60*x) + sin(70*x)

win = sin(linspace(0,pi,N));

Y = fftshift(abs(fft(fx)));

plot(1:N,Y)
%}

for ii = 1:30
    fprintf('| %1.3f  |  %1.3f  |\n', (ii-1)/ii,1/ii)
end