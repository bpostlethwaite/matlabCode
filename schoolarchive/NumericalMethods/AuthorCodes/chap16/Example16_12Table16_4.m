% Table 16.4 Example 16.12 : Adams-Moulton for a simple example
% Use dumb fixed point iteration to handle nonlinearity

clear all
format short e
h = [.2,.1,.05,.02,.01,.005,.002]; nh = length(h);

y0 = 1; t0 = 1; tf = 10; ye = 1/tf;
iters = 20;

% AM(1,1): backward Euler
y(1) = y0; t = t0;
for j = 1:nh
  k = h(j); N = (tf-t0) / k;
  for i = 1:N
    y(i+1) = y(i); % primitive initial guess  
    for ll=1:iters, y(i+1) = y(i) - k*y(i+1)^2; end
  end
  errAM1(j) = abs(ye - y(N+1));
end
rateAM1 = log2(errAM1(2:nh)./errAM1(1:nh-1)) ./ log2(h(2:nh)./h(1:nh-1));
errAM1
rateAM1

% AM(1,2) : implicit trap
y(1) = y0;
for j = 1:nh
  k = h(j); N = (tf-t0) / k;
  fip1 = -y(1)^2;
  for i = 1:N 
    fi = fip1; y(i+1) = y(i);
    for ll=1:iters, fip1 = -y(i+1)^2; y(i+1) = y(i) + .5*k*(fip1+fi); end
  end
  errAM2(j) = abs(ye - y(N+1));
end
rateAM2 = log2(errAM2(2:nh)./errAM2(1:nh-1)) ./ log2(h(2:nh)./h(1:nh-1));
errAM2
rateAM2

% AM(3,4)
y(1) = y0; f = zeros(1,4);
for j = 1:nh 
  k = h(j); N = (tf-t0) / k;
  % extra initial values
  for l=1:3, f(l+1) = -y(l)^2; y(l+1) = 1 / (t0 + k*l); end
  for i = 3:N  
    for l=1:3, f(l) = f(l+1); end   % advance the saved f values
    y(i+1) = y(i);
    for ll=1:iters
      f(4) = -y(i+1)^2; 
      y(i+1) = y(i) + k/24*(9*f(4)+19*f(3)-5*f(2)+f(1));
    end
  end
  errAM4(j) = abs(ye - y(N+1));
end
rateAM4 = log2(errAM4(2:nh)./errAM4(1:nh-1)) ./ log2(h(2:nh)./h(1:nh-1));
errAM4
rateAM4