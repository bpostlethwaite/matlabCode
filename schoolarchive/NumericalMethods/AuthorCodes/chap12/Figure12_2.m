% Figure 12.2 : Legendre polynomials

clf
x = -1:.01:1;
p0 = ones(size(x));
p1 = x;
p2 = .5*(3*x.^2 - 1);
p3 = 5/3*x.*p2 - 2/3*p1;
p4 = 7/4*x.*p3 - 3/4*p2;

plot(x,p0,x,p1,x,p2,x,p3,x,p4)
hold on
xlabel('x')
ylabel('\phi_j')
axis([-1 1 -1.1 1.1])