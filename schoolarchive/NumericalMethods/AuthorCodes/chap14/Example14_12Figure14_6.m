% Example 14.13 -- Figure 14.6 : 
% compare differentiation qualities of cheb and fft
clear all
close all

count = 0;
range = 10:10:100;
for n = range

%Chebyshev    
meth = 2;
count = count + 1;
a = -2*pi; b = 2*pi;
[xc, yc, ycp, De] = chebdif (@func,a,b,n,meth);
D2 = De*De; ycpp = D2*yc;

%err1 = abs(ycp - D*yc)
err_chebe(count) = max(abs(ycp - funcp(xc)));
err_chebe2(count) = max(abs(ycpp - funcpp(xc)));

%FFT
l = n/2; h = (b-a)/n;
x = a:h:b-h; x = x(:);
%x = a+h:h:b; x = x(:);
f = func(x);
f_hat = fft(f); % f_hat = f_hat(:);

r = 1; m2 = 0;
g_hat = (1i)^r * ([0:l-1 m2 -l+1:-1].^r)' .* f_hat;
g_hat = (2*pi/(b-a))^r * g_hat;
g = real(ifft(g_hat));
err_fft1(count) = max(abs(g - funcp(x)));

r = 2; m2 = l;
g_hat = (1i)^r * ([0:l-1 m2 -l+1:-1].^r)' .* f_hat;
g_hat = (2*pi/(b-a))^r*g_hat;
g = real(ifft(g_hat));
err_fft2(count) = max(abs(g - funcpp(x)));


end

figure(1)
subplot(2,1,1)
semilogy(range,err_chebe)
xlabel('n')
ylabel('cheb derivative error')
subplot(2,1,2)
semilogy(range,err_chebe2)
xlabel('n')
ylabel('2nd derivative error')

figure(2)
subplot(2,1,1)
semilogy(range,err_fft1)
xlabel('n')
ylabel('fft derivative error')
subplot(2,1,2)
semilogy(range,err_fft2)
xlabel('n')
ylabel('2nd derivative error')

figure(3)
subplot(2,2,1)
semilogy(range,err_fft1)
xlabel('n')
ylabel('fft 1st der error')
subplot(2,2,3)
semilogy(range,err_fft2)
xlabel('n')
ylabel('fft 2nd der error')
subplot(2,2,2)
semilogy(range,err_chebe)
xlabel('n')
ylabel('cheb 1st der error')
subplot(2,2,4)
semilogy(range,err_chebe2)
xlabel('n')
ylabel('cheb 2nd der error')