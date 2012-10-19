% Example 13.11 : checking DCT and DFT compression
close all
clear all
thr = .01; isw = 3;

x = 0:.001:1;
xx = 0:.0001:1;

if isw == 1
  y = exp(3*x) .* sin(100*pi*x.^2) ./ (1 + 20*x.^2);
  yy = exp(3*xx) .* sin(100*pi*xx.^2) ./ (1 + 20*xx.^2);
elseif isw == 2
  y = exp(-(x-.5).^2).* sin(20*pi*x) ./   (1 + 20*x.^2);
  yy = exp(-(xx-.5).^2).* sin(20*pi*xx) ./   (1 + 20*xx.^2);
elseif isw == 3
    y = log(1 + 2*pi*x);
    yy = log(1 + 2*pi*xx);    
end
subplot(2,2,1)
plot (xx,yy)

subplot(2,2,2)
plot (x,y)

% compress with fft

yhf = fft(y);
m = length(y); mf = m;
nyhf = norm(yhf)/ sqrt(m);
for i = 1:m
    if(abs(yhf(i)) < thr*nyhf), yhf(i) = 0; mf = mf - 1; end
end
ynf = ifft(yhf);
% number of remaining nonzero coefficients and compression error
mf  
errf = norm(real(ynf)-y)/norm(y)

subplot(2,2,3)
plot (x,real(ynf))

% compress with dct

yhc = dct(y);
m = length(y); mc = m;
nyhc = norm(yhc)/ sqrt(m);
for i = 1:m
    if(abs(yhc(i)) < thr*nyhc), yhc(i) = 0; mc = mc - 1; end
end
ync = idct(yhc);
% number of remaining nonzero coefficients and compression error
mc
errc = norm(ync-y)/norm(y)

subplot(2,2,4)
plot (x,ync)

figure(2)
plot(x,y,x,ync)
figure(1)





