% Example 13.8 -- Figure 13.6 : exp(ix) and roots of unity
clear all
theta = 0:.01:2*pi;
eix = exp(i*theta);
roots = exp(i*2*pi/8*[0:1:7]);

figure(1)
clf
plot(real(eix),imag(eix),real(roots),imag(roots),'rd');
hold on
xax = [-1.1 1.1]; yax = [0 0];
plot (xax,yax,'k',yax,xax,'k')
xlabel('\Re')
ylabel('\Im')


axis equal
axis ([-1.1 1.1 -1.1 1.1])


