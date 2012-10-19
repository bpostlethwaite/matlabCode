% Question 4
clear all
close all

t = [0;1;2];
tt = [0:0.1:2];
z = exp([0.1;0.9;2]);

b = log(z);
A = [ones(3,1),t];
u = A\b;
uu = [exp(u(1));u(2)];

b = z;
g = @(x,t)(x(1)*exp(x(2)*t));
J = @(x,t)([exp(x(2)*t), x(1)*t.*exp(x(2)*t)]);

x = uu; % Initial Guess from previous computations
x = [9;1];
iter = 0;

while true
    p = (J(x,t)'*J(x,t))\(J(x,t)'*(b - g(x,t)));
    x = x + p;
    iter = iter + 1;
    if norm(p) < 0.0001
       break
    end   
end

r1 = norm(g(uu,t) - b);
r2 = norm(g(x,t) - b);

fprintf(['linear inf norm is %1.6f \n'...
        'Gauss Norm is %1.6f \n'],r1,r2)

plot(t,z,'r*',tt,g(uu,tt),'s:',tt,g(x,tt),'--^')
xlim([-0.1,2.1])
legend('data','transformation fit','Gauss Newton','Location','Best')
title('Non-Linear Least Squares fit for Transformed equation vs Gauss Newton')
xlabel('t axis')
ylabel('data values')
text(0,7,sprintf('Inf norm of Linearized transform method: %1.4f',r1'),...
    'BackgroundColor',[.7 .9 .7])
text(0,6,sprintf('Inf norm of Gauss Newton method: %1.4f',r2'),...
    'BackgroundColor',[.7 .9 .7])

%figure(2)
%plot(t,g(uu,t) - b,t,g(x,t) - b)

sum(g(uu,t) - b)
sum(g(x,t) - b)