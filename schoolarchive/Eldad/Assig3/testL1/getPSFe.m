function[PSF] = getPSFe(n,sig)
% [PSF] = getPSFe(n,sig)
%

[x,y] = ndgrid(linspace(-1,1,n));

PSF = exp(-x.^2/sig - y.^2/sig);

