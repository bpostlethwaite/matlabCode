% Example 3.4 -- Figure 3.2

x = 0:.01:1.4;
y1 = x;
y2 = exp(-x);
plot (x,y1,x,y2)
hold on

n=8; xx = zeros(1,n); zer = xx;
xx(1) = 1;
for j=1:n-1
  xx(j+1) = exp(-xx(j));
end
itx(1:2:2*n-1) = xx;
itx(2:2:2*(n-1)) = xx(1:n-1);
ity(1) = 0;
ity(2:2:2*(n-1)) = xx(2:n);
ity(3:2:2*n-1) = xx(2:n);
plot(itx,ity,'r--')
hold on

plot(xx,zer,'k+')
hold on

xlabel('x')
gtext('y=x')
gtext('g(x)=e^x')
gtext('x_0')
gtext('x_1')
gtext('x_2')

