% Example 11.6 -- Figure 11.5 : constructing natural spline with one 
%   interior break point

x = [0,1,2];
y = [.1,.9,2];
y = [1.1,.9,2];

% initialize arrays
a = zeros(1,4); b = a;
% the simple coefficients
a(1) = y(1); b(1) = y(2);
% a 6 by 6 system for the rest
A = [1,1,1,0,0,0; 0,0,0,1,1,1; 1,2,3,-1,0,0; 0,2,6,0,-2,0; ...
     0,2,0,0,0,0; 0,0,0,0,2,6];
rhs = [y(2)-a(1);y(3)-b(1); 0; 0; 0; 0];
coefs = A \ rhs;
a(2) = coefs(1); a(3) = coefs(2); a(4) = coefs(3);
b(2) = coefs(4); b(3) = coefs(5); b(4) = coefs(6);

a
b
% plot
xx1 = 0:.01:1-.01;
xx2 = 1:.01:2;
xx = [xx1,xx2];
v1 = a(1) + a(2)*xx1 + a(3)*xx1.^2 + a(4)*xx1.^3;
v2 = b(1) + b(2)*(xx2-1) + b(3)*(xx2-1).^2 + b(4)*(xx2-1).^3;
v = [v1,v2];

plot(xx,v)
hold on
xlabel('x')
ylabel('v')

xb = 1; yb = .9;
plot(xb,yb,'ro')

gtext('s_0')
gtext('s_1')
