% Example 13.1 -- Figure 13.1 : somple low pass filtering

x = -pi:.001:pi;

figure(1)
clf

f = cos(3*x) - .5*sin(5*x) + .05*cos(104*x);
subplot(2,1,1)
plot(x,f)
xlabel('x')
ylabel('signal with noise')

v = cos(3*x) - .5*sin(5*x);
subplot(2,1,2)
plot(x,v)
xlabel('x')
ylabel('signal filtered')