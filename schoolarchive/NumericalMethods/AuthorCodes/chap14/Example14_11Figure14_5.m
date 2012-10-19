% Example 14.11 -- Figure 14.5 numerical differentiation at chebn points
clear all
close all

a = 0; b = pi; meth = 2;
count = 0;
for n = 5:5:40
  count = count + 1;
  [xc, yc, ycp, D] = chebdif (@func,a,b,n,meth);
  D2 = D*D; ycpp = D2*yc;
  err_cheb(count) = max(abs(ycp - funcp(xc)));
  err_cheb2(count) = max(abs(ycpp - funcpp(xc)));
end

subplot(2,1,1)
semilogy(5:5:40,err_cheb)
xlabel('n')
ylabel('derivative error')
subplot(2,1,2)
semilogy(5:5:40,err_cheb2)
xlabel('n')
ylabel('2nd derivative error')