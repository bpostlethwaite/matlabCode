function [RdPP,RdSP,RdPS,RdSS,TuPP,TuSP,TuPS,TuSS]=rtcoef(a1,b1,r1,a2,b2,r2,p);

% FUNCTION [RdPP,RdSP,RdPS,RdSS,TuPP,TuSP,TuPS,TuSS] = rtcoef(a1,b1,r1,a2,b2,r2,p);
% Formulae for reflection (downward incidence) and transmission 
% coefficients (upward incidence) from Aki and Richards 1980, 
% pp 149-151. 

p2=p.*p;
qa1=sqrt(1/(a1*a1)-p2);
qa2=sqrt(1/(a2*a2)-p2);
qb1=sqrt(1/(b1*b1)-p2);
qb2=sqrt(1/(b2*b2)-p2);

a=r2*(1-2*b2*b2*p2)-r1*(1-2*b1*b1*p2);
b=r2*(1-2*b2*b2*p2)+2*r1*b1*b1*p2;
c=r1*(1-2*b1*b1*p2)+2*r2*b2*b2*p2;
d=2*(r2*b2*b2-r1*b1*b1);

E=b.*qa1+c.*qa2;
F=b.*qb1+c.*qb2;
G=a-d.*qa1.*qb2;
H=a-d.*qa2.*qb1;
D=E.*F+G.*H.*p2;

RdPP=((b.*qa1-c.*qa2).*F-(a+d.*qa1.*qb2).*H.*p2)./D;
RdSP=-2*qa1.*(a.*b+c.*d.*qa2.*qb2).*p.*a1./(b1*D);
RdPS=-2*qb1.*(a.*b+c.*d.*qa2.*qb2).*p.*b1./(a1*D);
RdSS=-((b.*qb1-c.*qb2).*E-(a+d.*qa2.*qb1).*G.*p2)./D;

TuPP=2*r2*qa2.*F.*a2./(a1*D);
TuSP=-2*r2*qa2.*G.*p*a2./(b1*D);
TuPS=2*r2*qb2.*H.*p*b2./(a1*D);
TuSS=2*r2*qb2.*E.*b2./(b1*D);
