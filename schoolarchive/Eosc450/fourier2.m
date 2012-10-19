clear all
close all
x=-pi:0.01:pi;
m=[2,5,10,100];
ax=[-pi,pi,0,5];

for ii = 1:4
 G = cos(x);
    for n = 2:m(ii)  
        f= (cos((2*n-1).*x) )/(2*n-1)^2;
        G = G + f;
  
    end
 F = pi/2 -(4/pi)*G;

 subplot(2,2,ii);
 plot(x,F);
 axis(ax);
 title(['value of n is ',num2str(m(ii))])
 hold on
 plot(x,abs(x),'--r');
 
end

