% Figure 10.3
%
% plot ith Lagrange polynomial of degree n based on
% n+1 data abscissae xx
clf

xx = [1,2,4];
n = 2;
a = .2;
b = 5.2;
t = a:.01 : b;

for i = 0:n
  xi = xx(i+1);
  denom = 1;
  for k = 0:n
    if k ~= i
      denom = denom * (xi - xx(k+1));
    end
  end

  for l = 1:length(t)
    li(i+1,l) = 1/denom;
    for k = 0:n
      if k ~= i
        li(i+1,l) = li(i+1,l) * (t(l) - xx(k+1));
      end
    end
  end
end

plot (t,li)
hold on
xlabel('x')
ylabel('Quadratic Lagrange polynomials')
gtext('L_0(x)')
gtext('L_1(x)')
gtext('L_2(x)')

