% Question 5
clear all
close all

n = 33;

t = linspace(0.9,4.1,n);
tt = linspace(0.9,4.1,n*10);

q = -11 + 55/3*t - 17/2*t.^2 + 7/6*t.^3;

b = q.*(1 + 0.3*randn(1,n));



pfcoef = polyfit(t,b,22);
PF = polyval(pfcoef,t);
SP = spline(t,b,tt);
lscoef = lsfit(t,b, 4);
A = [ones(length(tt),1),tt(:),tt(:).^2,tt(:).^3];
LS = A*lscoef;


figure(1)
subplot(3,1,1)
    plot(t,b,'*',t,PF)
    title('Polyfit of degree 32')
    legend('Data Points','polyfit interpolation','Location','Best')
    ylim([-4,4])
subplot(3,1,2)
    plot(t,b,'*',tt,SP)
    title('Spline interpolation')
    legend('Data Points','Spline interpolation','Location','Best')
subplot(3,1,3)
    plot(t,b,'*',tt,LS)
    title('lsfit with cubic polynomial')
    legend('Data Points','lsfit cubic interpolation','Location','Best')