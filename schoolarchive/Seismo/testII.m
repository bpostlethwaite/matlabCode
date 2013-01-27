clear all
close all

w = 0.5;
t = 0:0.01:2;
A = 0.1;


u = A*sin(2*pi*w*t);

v = 2*pi*w*A*cos(2*pi*w*t);

a = -4*(pi^2)*(w^2)*A*sin(2*pi*w*t);

plot(t,u,t,v,t,a)
legend('Displacement','Velocity','Acceleration','Location','Best')