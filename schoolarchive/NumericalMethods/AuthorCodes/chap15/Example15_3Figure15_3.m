% Example 15.3 -- Figure 15.3 : quadrature errors for int_0^1 exp(-x^2)
clear all
close all

exact = 0;
%trap approximations
count = 0;
for r = 4:2:16
    h = 1/r; x = 0:h:1; count = count + 1;
    f = exp(-x.^2);
    itr(count) = h*sum(f) - h/2*(f(1)+f(end));
    errt(count) = abs(exact - itr(count));
end

%Simpson approximations
count = 0;
for r = 4:2:32
    h = 1/r; x = 0:h:1; count = count + 1;
    f = exp(-x.^2);
    isp(count) = h/3*(f(1)+f(end));
    isp(count) = isp(count)+ 2*h/3* sum(f(3:2:r-1)) + 4*h/3* sum(f(2:2:r));
    errs(count) = abs(exact - isp(count));
end

exact = errs(end);
ert = abs(exact - itr);
ers = abs(exact - errs(1:7));

semilogy (4:2:16,ert,4:2:16,ers)
hold on
legend('error_{trap}','error_{Simp}')
xlabel('r')

