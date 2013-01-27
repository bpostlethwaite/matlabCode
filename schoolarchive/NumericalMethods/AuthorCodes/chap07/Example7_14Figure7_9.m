% Example 7.14 -- Figure 7.9 : polynomials of CG iteration

% small 3 x 3 problem
A=[7 3 1; 3 10 2; 1 2 15];
b=[28;31;22];

eigA=eig(A);
r0=b; p0=b;
delta0=norm(r0)^2;
s0=A*p0;
alpha0=delta0/(p0'*s0);
x1=zeros(3,1)+alpha0*p0;
r1=b-alpha0*s0;
delta1=norm(r1)^2;
p1=r1+delta1/delta0*p0;
s1=A*p1;
alpha1=delta1/(p1'*s1);
x2=x1+alpha1*p1;

r2=r1-alpha1*s1;
delta2=norm(r2)^2;
p2=r2+delta2/delta1*p1;
s2=A*p2;
alpha2=delta2/(p2'*s2)
x3=x2+alpha2*p2


x=0:.01:20;

% p1:
c10=1;
c11=-alpha0;
p1=c10+c11*x;

% p2:
c20 = 1                                      % const (always 1)
c21 = -(alpha0+alpha1+alpha1*delta1/delta0)  % coef of x
c22 = alpha0*alpha1                          % coef of x^2
p2=c20+c21*x+c22*x.^2;

% p3:
c30 = 1;
c31 = -alpha0-alpha1*(1+delta1/delta0)-alpha2*(1+delta2/delta0+delta2/delta1);
c32 = alpha0*alpha1+alpha2*(alpha0+alpha1+alpha1*delta1/delta0)+alpha0*alpha2*delta2/delta1;
c33=-alpha0*alpha1*alpha2    % should be equal to: -1/(-eigA(1)*eigA(2)*eigA(3));
p3=c30+c31*x+c32*x.^2+c33*x.^3;

% Can also do this and then we are sure we're OK...
p33=(x-eigA(1)).*(x-eigA(2)).*(x-eigA(3))/(-eigA(1)*eigA(2)*eigA(3));

xx=eigA;
p1_eig=c10+c11*xx;
p2_eig=c20+c21*xx+c22*xx.^2;
p3_eig=c30+c31*xx+c32*xx.^2+c33*xx.^3;

plot(x,p1);
hold on
plot(x,p2);
plot(x,p3)
plot(xx,p1_eig,'o')
plot(xx,p2_eig,'+')
plot(xx,p3_eig,'*')

line([min(x),max(x)],[0 0]);
