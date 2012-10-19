function tol = LCurveTolPlot(A, b)
%
%            tol = LCurveTolPlot(A, b);
%
%  Given a psfMatrix A, and RHS image b, this function uses
%  an L-Curve to help find a tolerance for the preconditioner.
%

%  J. Nagy  12/30/01

P = psfPrec(A, b, 0);
e = 1./P.matdata;
e = reshape(e', prod(size(e)), 1);
[e, idx] = sort(abs(e));
e = flipud(e);

fb = reshape(fft2(b)', prod(size(b)), 1);
fb = abs(flipud(fb(idx))) / sqrt(prod(size(b)));

eta = fb .* fb;
x = eta ./ (e.^2);

rho = zeros(size(x));
sigma = rho;
rho(1) = sum(eta(2:end));
sigma(1) = x(1);

for k = 2:length(x)-1
  rho(k) = rho(k-1) - eta(k);
  sigma(k) = sigma(k-1) + x(k);
end

figure
h = gcf;
loglog(rho, sigma)
xlabel('norm of residual')
ylabel('norm of solution')
hold on
drawnow

disp('click near the corner to zoom in, and hit the space bar to continue')
zoom
pause

disp(' ')
disp('Use the mouse to locate corner of L-curve')
[rr, xx] = ginput(1);
str{1} = ' ';

k = round( mean(sum(rho > rr), sum(sigma < xx)) );
tol = e(k);
str{2} = sprintf('tol = %f', tol);

plot(rho(k), sigma(k), 'o')
legend(str)

NewPoint = input('Choose a different corner? [n] ', 's');
if isempty(NewPoint)
  NewPoint = 'n';
end
l = 2;
pt = ['go'; 'ro'; 'co'; 'mo'; 'yo'; 'ko'];
ptl = 1;
while NewPoint == 'y'
  l = l + 1;
  figure(h)
  disp('   Use the mouse to locate corner of L-curve')
  [rr, xx] = ginput(1);
  k = round( mean(sum(rho > rr), sum(sigma < xx)) );
  tol = e(k);
  str{l} = sprintf('tol = %f', tol);
  plot(rho(k), sigma(k), pt(ptl,:))
  legend(str)
  
  NewPoint = input('Choose a different corner? [n] ', 's')
  if isempty(NewPoint)
    NewPoint = 'n';
  end
  if ptl == 6
    ptl = 1;
  else
    ptl = ptl + 1;
  end
end
