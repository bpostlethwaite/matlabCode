% examples of Chebyshev interpolation

clear all
which = 2; % should be same value in func

if which == 1  % Runge example
    a = -1; b = 1; 
    nn = 10:10:170;
elseif which == 2 % highly oscillatory
    a = 0; b = 1;
    nn = 50:10:200;
elseif which == 3 % hat function
    a = 0; b = 2*pi;
    nn = 1:2:31; 
elseif which == 4 % example 13.3
    a = -1; b = 2;
    nn = 3:2:31; 
elseif which == 5 % example 13.4
    a = -1; b = 2;
    nn = [3,7,15,31,63,127]; 
elseif which == 6 % example 13.1
    a = -pi; b = pi;
    nn = 10:10:200; nn = nn - 1; 
end
xe = a:.001:b;

for j = 1:length(nn)
  n = nn(j);  
  [xc,yc,rho] = chebl(@func,a,b,n);
  
  fe = func(xe);
  pe = lageval(xe,xc,yc,rho);
  err(j) = max(abs(fe-pe));
  err_max(j) = max(abs(fe-pe) ./ (abs(fe)+1));
end

figure(1)
clf
plot(xe,fe)
xlabel('x')
ylabel('f')
pe = lageval(xe,xc,yc,rho);
figure(2)
clf
plot(xe,pe)
xlabel('x')
ylabel('p_{n}')
    
figure(3)
clf
semilogy(nn,err)
xlabel('n')
ylabel('max error')

if which == 2
figure(4)
clf
subplot(2,1,1)
plot(xe,pe)
xlabel('x')
ylabel('p_{200}')
subplot(2,1,2)
semilogy(nn,err)
xlabel('n')
ylabel('max error')
end