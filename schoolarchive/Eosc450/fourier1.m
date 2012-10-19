clear all
close all
x=-pi:0.01:pi;
m=[2,5,10,100];
ax=[-pi,+pi,-2,12];

for ii = 1:4
G = -(0.5).*(cos(x)-sin(x));  
    for n = 2:m(ii)  
        f= (((-1)^n)/(n^2+1)).*(cos(n*x)-n*sin(n*x));  
        G = G + f;
  
    end
F = sinh(pi)/pi + (2*sinh(pi)/pi)*G;

 subplot(2,2,ii);
 plot(x,F);
 axis(ax);
 title(['value of n is ',num2str(m(ii))])
 hold on
 plot(x,exp(x),'--r');
end