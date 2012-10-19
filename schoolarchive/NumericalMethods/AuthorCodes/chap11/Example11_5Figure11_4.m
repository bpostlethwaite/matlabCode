% Example 11.5 -- Figure 11.4 : spline interpolation of Runge example

clear all
clf
h = 2/19;
xx = -1:h:1;
yy=1 ./(1+25*xx.^2);

x = -1:.001:1;
tt = .5*(xx+1);
t = .5*(x+1);

u = 1 ./(1+25*x.^2);
y = spline(tt,yy,t);

plot(tt,yy,'o',t,u,t,y)
hold on
xlabel('t = .5*(x+1)')
ylabel('v')
legend('interp points','exact','spline') 
axis([0 1 0 1.05])

max(abs(u-y))