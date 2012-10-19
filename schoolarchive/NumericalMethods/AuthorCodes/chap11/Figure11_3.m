% Figure 11.3 : matching two cubics with different 3rd deriv.

clear all
close all
A = [.05, .05^2, .05^3; .1, .1^2, .1^3; -.05, (-.05)^2, 7*(-.05)^3];
rhs = [1.5-1; 1-1; .25-1];

x = A\rhs;

ai = [1;x];
aim = ai; aim(4) = 7*ai(4);

x1 = -.1:.01:0;
v1 = aim(1) + aim(2)*x1 + aim(3)*x1.^2 + aim(4)*x1.^3;
x2 = 0:.01:.14;
v2 = ai(1) + ai(2)*x2 + ai(3)*x2.^2 + ai(4)*x2.^3;
clf
plot(x1,v1,'b',x2,v2,'g',0,1,'rd')
hold on
xlabel('x')
ylabel('v')
gtext('s_i')
gtext('s_{i-1}')
gtext('x_i')