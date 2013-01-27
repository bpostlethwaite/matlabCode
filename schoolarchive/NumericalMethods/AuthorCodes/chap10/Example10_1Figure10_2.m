% Example 10.1 -- Figure 10.2

clf
x1 = [1,1.5,2,3,4]; x2 = x1;
c1 = [2,-1]; 
c2 = [-2,12,-7]/3;
y1 = c1(1)*x1 + c1(2);
y2 = c2(1)*x2 .* x2 + c2(2)*x2 + c2(3);
xx1 = .3:.1:4.3;
zz1 = c1(1)*xx1 + c1(2);
xx2 = .3:.1:5.5;
zz2 = c2(1)*xx2 .* xx2 + c2(2)*xx2 + c2(3);

plot(x1,y1,'md',x2,y2,'ro',xx1,zz1,'b',xx2,zz2,'g')
hold on
xlabel('x')
ylabel('p')
gtext('p_1(x) = 2x-1')
gtext('p_2(x) = (-2x^2+12x-7)/3')