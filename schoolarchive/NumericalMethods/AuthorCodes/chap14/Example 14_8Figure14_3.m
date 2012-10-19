% Example 14.6 -- Figure 14.3 : differentiating noisy data

x = 0:.01:2*pi;
l = length(x);
sinx = sin(x);
sinp = (1+.01*randn(1,l)).*sinx;

cosx = (sinx(3:l)-sinx(1:l-2))/.02;
cosp = (sinp(3:l)-sinp(1:l-2))/.02;
err_f = max(abs(sinx-sinp))
err_fp = max(abs(cosx-cosp))

subplot(1,2,1)
plot(x,sinp,x,sinx,'r')
xlabel('x')
%title('sin (x) with 1% noise')

subplot(1,2,2)
plot(x(2:l-1),cosp,x(2:l-1),cosx,'r')
xlabel('x')
%title('cos (x) by noisy numerical differentiation')