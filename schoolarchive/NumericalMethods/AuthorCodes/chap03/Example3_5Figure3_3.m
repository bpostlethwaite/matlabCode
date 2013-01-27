% Example 3.5 -- Figure 3.3

x = [0 10];
plot(x,x,'g')
hold on
axis ([0 10 0 14])

x = 0:.01:10;
y = 2* cosh(x/4);
plot(x,y,'b')
xlabel('x')

% find the two points where these curves meet and mark them
a = 2;
for k = 1:4
    a = a - (2*cosh(a/4) - a)/(.5*sinh(a/4)-1);
end
b = 8;
for k = 1:4
    b = b - (2*cosh(b/4) - b)/(.5*sinh(b/4)-1);
end
plot([a b],[a b],'g')

gtext('y=x')
gtext('y=2cosh(x/4)')
