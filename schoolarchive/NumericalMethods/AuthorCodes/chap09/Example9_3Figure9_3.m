% Example 9.3 -- Figure 9.3
% Program to solve  u'' + e^u = 0, u(0)= u(1) = 0

J = 24;  % How many interior mesh points
h = 1 / (J+1);
t = 0:h:1;  % spatial mesh

% initial guess
alfa = input ('enter alfa between 0 and 100 : ');
u0e = alfa * t.*(1-t);
u0 = u0e(2:J+1);

% find and plot corresponding solution
[u,k] = newtons('funcode',u0,1.e-8,20);
figure(1)
plot(t,[0;u;0]')
hold on
