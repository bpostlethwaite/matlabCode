% Table 16.3 Example 16.12 : Adams-Bashforth for a simple example

clear all
format short e
h = [.2,.1,.05,.02,.01,.005,.002]; nh = length(h);

y0 = 1; t0 = 1; tf = 10; ye = 1/tf;

% AB(1,1): forward Euler
y(1) = y0; t = t0;
for j = 1:nh
  k = h(j); N = (tf-t0) / k;
  for i = 1:N, y(i+1) = y(i) - k*y(i)^2; end
  errAB1(j) = abs(ye - y(N+1));
end
rateAB1 = log2(errAB1(2:nh)./errAB1(1:nh-1)) ./ log2(h(2:nh)./h(1:nh-1));
errAB1
rateAB1

% AB(2,2)
y(1) = y0;
for j = 1:nh
  k = h(j); N = (tf-t0) / k;
  fi = -y(1)^2; y(2) = 1/(t0+k);         % extra initial value
  for i = 2:N 
    fim1 = fi; fi = -y(i)^2; y(i+1) = y(i) + .5*k*(3*fi-fim1); 
  end
  errAB2(j) = abs(ye - y(N+1));
end
rateAB2 = log2(errAB2(2:nh)./errAB2(1:nh-1)) ./ log2(h(2:nh)./h(1:nh-1));
errAB2
rateAB2

% AB(4,4)
y(1) = y0; f = zeros(1,4);
for j = 1:nh 
  k = h(j); N = (tf-t0) / k;
  % extra initial values
  for l=1:3, f(l+1) = -y(l)^2; y(l+1) = 1 / (t0 + k*l); end
  for i = 4:N  
    for l=1:3, f(l) = f(l+1); end  
    f(4) = -y(i)^2; 
    y(i+1) = y(i) + k/24*(55*f(4)-59*f(3)+37*f(2)-9*f(1));
  end
  errAB4(j) = abs(ye - y(N+1));
end
rateAB4 = log2(errAB4(2:nh)./errAB4(1:nh-1)) ./ log2(h(2:nh)./h(1:nh-1));
errAB4
rateAB4