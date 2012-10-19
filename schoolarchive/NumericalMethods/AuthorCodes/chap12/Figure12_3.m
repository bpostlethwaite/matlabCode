% Figure 12.3 : Chebyshev polynomials

clf
x = -1:.01:1;
p0 = ones(size(x));
p1 = x;
p2 = 2*x.*p1 - p0;
p3 = 2*x.*p2 - p1;
p4 = 2*x.*p3 - p2;
p5 = 2*x.*p4 - p3;
p6 = 2*x.*p5 - p4;

plot(x,p4,x,p5,x,p6)
hold on
xlabel('x')
ylabel('T_j')
gtext('T_4')
gtext('T_5')
gtext('T_6')