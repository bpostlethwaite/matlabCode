% Example 15.4 -- Figure 15.4 : quadrature error for int_10^10 exp(-x^2)
%     using composite rules
clear all
figure(1)
clf

exact = 0;
%trap approximations
count = 0;
for r = 32:8:128
    h = 20/r; x = -10:h:10; count = count + 1;
    f = exp(-x.^2);
    itr(count) = h*sum(f) - h/2*(f(1)+f(end));
    errt(count) = abs(exact - itr(count));
end

%Simpson approximations
count = 0;
for r = 32:8:64
    h = 20/r; x = -10:h:10; count = count + 1;
    f = exp(-x.^2);
    isp(count) = h/3*(f(1)+f(end));
    isp(count) = isp(count)+ 2*h/3* sum(f(3:2:r-1)) + 4*h/3* sum(f(2:2:r));
    errs(count) = abs(exact - isp(count));
end

exact = itr(end);
ert = abs(exact - itr(1:5));
ers = abs(exact - isp(1:5));

semilogy (32:8:64,ert,32:8:64,ers)
hold on
legend('error_{trap}','error_{Simp}')
xlabel('r')

