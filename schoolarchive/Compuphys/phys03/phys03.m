close all
clear all
load phys03.mat
x = [0.1 1 100];
for ii = 1:3

h = A{ii}(:,5);
forward = A{ii}(:,6);
center = A{ii}(:,7);
extrap = A{ii}(:,8);
alg1 = B{ii}(:,6);
alg2 = B{ii}(:,7);

E_prec1 = log10(abs(2*1e-16*(-sin(x(ii)))./A{ii}(:,1)));
E_prec2 = log10(abs(2*1e-16*(-cos(ii))./A{ii}(:,1).^2));

figure(1)
subplot(3,1,ii)
plot(h,forward,'b*',h,center,'r*',h,extrap,'k*');
hold on
plot(h,E_prec1,'g','LineWidth',3);
xlabel('Log h')
ylabel(' Log Error')
title(sprintf('LogLog plot of relative error vs h, for 1st derivatives &  x = %3.1f',x(ii)))
legend('Forward Difference','Centred Difference','Extrapolated Difference','Precision Error','Location','Best')

figure(2)
subplot(3,1,ii)
plot(h,alg1,'b*',h,alg2,'r+');
hold on
plot(h,E_prec2,'g','LineWidth',3);
xlabel('Log h')
ylabel(' Log Error')
title(sprintf('LogLog plot of relative error vs h, for 2nd derivatives & x = %3.1f',x(ii)))
legend('Algorithm 1','Algorithm 2','Precision Error','Location','Best')

end