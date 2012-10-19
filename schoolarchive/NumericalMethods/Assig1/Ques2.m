clear all
close all

n = 9;
x = 1.84:.001:2.16;
c = [-512,2304,-4608,5376,-4032,2016,-672,144,-18,1];
p = c(n+1);

for i = n:-1:1
        p = p.*x + c(i);
end

pd = (x-2).^9;
relErr = abs(pd-p)./p;
figure(1)
  subplot(2,1,1)
    plot(x,p)
    title('Nested polynomial calculation')
    xlim([1.94,2.08])
    legend('Nested','Location','NorthWest')
    xlabel('x')
    ylabel('f(x)')
  subplot(2,1,2)
    plot(x,pd)
    title('Direct polynomial calculation')
    legend('Direct','Location','NorthWest')
    xlim([1.94,2.08])
    xlabel('x')
    ylabel('f(x)')
    
figure(2)
    plot(x,relErr)
    title('relative error of Nested vs direct calculation')
    xlabel('x')
    ylabel('relative error')
    ylim([-3,3])
