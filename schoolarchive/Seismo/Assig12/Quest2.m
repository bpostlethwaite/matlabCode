clear all
close all
modes = 1;
beta1 = 3500;
beta2 = 4500;
c = linspace(3000,3400,1000)';
p = 1./c;
alpha = [1.5:0.01:2] * beta1;
beta = beta1;

for jj = 1:length(p)
f1(jj,:) = (2*p(jj)^2 - 1./beta^2).^2;
f2(jj,:) = 4*p(jj)^2*sqrt(p(jj)^2 - 1./alpha.^2).*sqrt(p(jj).^2 - 1./beta^2);
end

for jj = 1:length(alpha)
A(jj,:) = getequal(f1(:,1),f2(:,jj),p,modes);
end

C = 1./A;

Egg = C./beta;

figure(1)
plot(f1)
hold on
plot(f2,':')

figure(2)
plot(alpha./beta,Egg);