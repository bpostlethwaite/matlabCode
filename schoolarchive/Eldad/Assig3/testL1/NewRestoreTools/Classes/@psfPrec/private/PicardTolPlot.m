function tol = PicardTolPlot(A, b)
%
%             tol = PicardTolPlot(A, b);
%
%  Given a psfMatrix A, and a right hand side image b,
%  this function plots the Picard condition, and allows
%  the user to graphically choose a tolerance for use in
%  the preconditioner.
%

%  J. Nagy  12/30/01

P = psfPrec(A, b, 0);
e = 1./P.matdata;
e = abs(reshape(e', prod(size(e)), 1));
[e, idx] = sort(e);
e = flipud(e);

b = reshape( fftn(b)', prod(size(b)), 1 );
b = abs( flipud( b(idx) ) ) / sqrt(prod(size(b)));

figure 
h = gcf;
semilogy(e)
hold on
semilogy(b, 'r')
legend('Eigenvalues', 'Fourier coefficients')
drawnow

SmoothPlot = input('Smooth the Fourier coefficeints? [n] ', 's');
if isempty(SmoothPlot)
  SmoothPlot = 'n';
end
l = 1;
while SmoothPlot == 'y' 
  hold off
  figure(h)
  semilogy(e)
  hold on
  semilogy(medfilt1(b,l*10), 'r')
  legend('Eigenvalues', 'Fourier coefficients')
  drawnow
  if l > 5
    break
  end
  SmoothPlot = input('Smooth the Fourier coefficeints? [n] ', 's');
  if isempty(SmoothPlot)
    SmoothPlot = 'n';
  end
  l = l + 1;
end


disp(' ')
disp('Use the mouse to locate where Fourier coefficients level off ')

[k, y] = ginput(1);
k = round(k);
tol = e(k);

plot(k, e(k), 'o')
plot([k, k], [min(b), e(k)], 'k--')
text(k, e(k), sprintf('  tol = %f', tol))

NewPoint = input('Choose a different point? [n] ', 's');
if isempty(NewPoint)
  NewPoint = 'n';
end

while NewPoint == 'y'
  figure(h)
  disp('Use the mouse to locate where the Fourier coefficients level off')
  [k, y] = ginput(1);
  k = round(k);
  tol = e(k);
  plot(k, e(k), 'o')  
  plot([k,k], [min(b), e(k)], 'k--')
  text(k, e(k), sprintf('  tol = %f', tol))

  NewPoint = input('Choose a different point? [n] ', 's');
  if isempty(NewPoint)
    NewPoint = 'n';
  end
end


  
