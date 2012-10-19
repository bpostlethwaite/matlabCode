% Example 15.10 -- Figure 15.5 : discontinuous integrand figure 15.5

clear all
close all
figure(1)
x = 0:.01*pi:2*pi;
f = sin(x);
ii = find(x > pi/2);
f(ii:end) = cos(x(ii:end));
plot (x,f)
xlabel('x')
ylabel('f')
axis([-.5 7 -1.1 1.1])

figure(2)
x1 = 0:.01*pi:pi/2;
f1 = sin(x1);
x2 = [pi/2, pi/2, pi/2, pi/2, pi/2];
f2 = [1,.75, .5, .25, 0];
x3 = pi/2:.01*pi:2*pi;
f3 = cos(x3);
plot (x1,f1,'b',x2,f2,'r',x3,f3,'b')
xlabel('x')
ylabel('f')
axis([-.5 7 -1.1 1.1])