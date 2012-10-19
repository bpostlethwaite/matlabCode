% Example 5.22 -- Figure 5.8

clear all
d = .75;
A = [1,d;d,1];
t = 0:.01:10; m = length(t);
x(1,1:m) = sin(t);
x(2,1:m) = cos(t);
y = A*x;
plot (y(1,:),y(2,:),'b', x(1,:),x(2,:),'r:')
xlabel('x_1')
ylabel('x_2')