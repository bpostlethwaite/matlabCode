% Figure 6.3 : least squares data fitting by a straight line

t = 1:1:25;
b = [1.0,1.1,1.0,1.2,1.2,1.3,1.4,1.5,1.3,1.4,1.6,1.9,1.6,1.8,1.8,1.9,...
    2.0,2.0,2.1,2.2,2.3,2.2,2.4,2.2,2.5];

coefs = lsfit(t,b,2);

v = coefs(1) + coefs(2)*t;

figure(2)
plot(t,b,'go',t,v)
hold on
xlabel('t')
legend('data','v')