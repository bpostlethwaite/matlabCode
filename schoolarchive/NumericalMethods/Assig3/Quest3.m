% Question 3
clear all
close all

t = [0;1;2];
tt = [0:0.1:2];
z = exp([0.1;0.9;2]);
b = log(z);
A = [ones(3,1),t];
u = A\b;

uu = exp(u);
ut = uu(1)*exp(u(2)*tt);
plot(t,z,'r*',tt,ut)
xlim([-0.1,2.1])
legend('data','fit','Location','Best')