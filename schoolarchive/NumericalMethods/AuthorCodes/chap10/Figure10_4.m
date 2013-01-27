% Figure 10.4
% plot ith Lagrange polynomial of degree n based on n+1 data abscissae xx
clear all
clf
xx = 1:1:6;
i = 2;
n = length(xx) - 1;
xi = xx(i+1);
xmin = min(xx); xmax = max(xx);
a = xmin - 0.01;
b = xmax + 0.01;
t = a:.01 : b;

denom = 1;
for k = 0:n
  if k ~= i
    denom = denom * (xi - xx(k+1));
  end
end

for l = 1:length(t)
  li(l) = 1/denom;
  for k = 0:n
    if k ~= i
      li(l) = li(l) * (t(l) - xx(k+1));
    end
  end
end

plot (t,li)
hold on
%title('Lagrange polynomial')
xlabel('x')
ylabel('L_2')
x1 = [.001,6.999];
y1 = [0,0];
plot(x1,y1,'k')
yy = zeros(size(xx));
plot (xx,yy,'r+')

