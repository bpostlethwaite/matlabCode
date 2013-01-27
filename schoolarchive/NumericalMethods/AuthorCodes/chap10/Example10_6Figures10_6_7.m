% Example 10.6 -- Figures 10.6 and 10.7 : interpolate the Runge example
clear all
close all

%quartic interpolant at equidistant points
xi4 = [-1,-.5,0,.5,1]; yi4 = 1 ./(1+25*xi4.^2);
coef4 = divdif(xi4,yi4);
%evaluate quartic
x = -1.:.002:1.;
ye = 1 ./ (1+25*x.^2); 
y4 = evalnewt(x,xi4,coef4);
figure(1)
plot (x,ye,'b',x,y4,'r',xi4,yi4,'co')
hold on
xlabel('x')
ylabel('p')
axis([-1.1,1.1,-0.5,1.1])

%9deg interpolant
h = 2/9;
xi9 = -1:h:1-10*eps; yi9 = 1 ./ (1+25*xi9.^2);
coef9 = divdif(xi9,yi9);
%evaluate 9deg poly
y9 = evalnewt(x,xi9,coef9);
plot (x,y9,'g',xi9,yi9,'cd')

%19deg interpolant
h = 2/19;
xi19 = -1:h:1-10*eps; yi19 = 1 ./ (1+25*xi19.^2);
coef19 = divdif(xi19,yi19);
%evaluate 9deg poly
y19 = evalnewt(x,xi19,coef19);
figure(2)
plot (x,ye,'b',x,y19,'m',xi19,yi19,'co')
hold on
xlabel('x')
ylabel('p')
axis([-1.1,1.1,-1,9])

%quartic interpolant at Chebyshev points
i4 = 0:1:4; xi4 = cos(pi/(2*(4+1))* (2*i4+1)); yi4 = 1 ./(1+25*xi4.^2);
coef4 = divdif(xi4,yi4);
%evaluate quartic
x = -1.:.002:1.;
ye = 1 ./ (1+25*x.^2); 
y4 = evalnewt(x,xi4,coef4);
figure(3)
plot (x,ye,'b',x,y4,'r',xi4,yi4,'co')
hold on
xlabel('x')
ylabel('p')
axis([-1.1,1.1,-0.5,1.1])

%9deg interpolant
i9 = 0:1:9;
xi9 = cos(pi/(2*(9+1))* (2*i9+1)); yi9 = 1 ./ (1+25*xi9.^2);
coef9 = divdif(xi9,yi9);
%evaluate 9deg poly
y9 = evalnewt(x,xi9,coef9);
plot (x,y9,'g',xi9,yi9,'cd')

%19deg interpolant
i19 = 0:1:19;
xi19 = cos(pi/(2*(19+1))* (2*i19+1)); yi19 = 1 ./ (1+25*xi19.^2);
coef9 = divdif(xi19,yi19);
%evaluate 9deg poly
y19 = evalnewt(x,xi19,coef9);
figure(4)
plot (x,ye,'b',x,y19,'m',xi19,yi19,'co')
hold on
xlabel('x')
ylabel('p')
axis([-1.1,1.1,-0.5,1.1])

