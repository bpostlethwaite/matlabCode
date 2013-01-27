% Example 9.1 -- Figure 9.1 : simple plots

x = -0.4:.01:1.1;
y = x.^2 -2*x + 1;

t = 0:.01:2*pi;
xx = sin(t);
yy = cos(t);

plot(x,y,xx,yy)
axis([-1.8 1.8 -1 2])
xlabel('x_1')
ylabel('x_2')