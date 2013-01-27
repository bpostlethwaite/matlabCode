% Example 10.9 -- Figure 10.10 : Hermite interpolation of log (x)
clear all
clf

% by monomial
format long g
c = hermite(1,2,0,1,log(2),.5)
x = .5:.01:2.5;
y1 = log(x);
y2 = c(4)*x.^3 + c(3)*x.^2 + c(2)*x.^1 + c(1)*ones(size(x));
plot(x,y1,'g',x,y2,'b')
hold on
xlabel('x')
ylabel('y')
legend('ln','Hermite cubic')
plot([1,2],[0,log(2)],'go')