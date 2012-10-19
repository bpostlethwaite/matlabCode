
m=[2,5,10,50];
x=linspace(-pi,pi);



for jj = 1:4
    
    Y =-0.5*(cos(x)-n*sin(x));
    
for n = 2:m(jj)
    
    G =( (((-1)^n)/(n^2+1))*(cos(n*x)-n*sin(n*x)));
    Y = Y + G;
    
   
F=sinh(pi)/pi+2*sinh(pi)/pi * Y;

subplot(2,2,jj);
  plot(x,F);
  title(['for n = ',num2str(m(jj)),]);


end

end
