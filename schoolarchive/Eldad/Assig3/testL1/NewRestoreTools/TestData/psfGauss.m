function PSF = psfGauss(dim, sigma)
%
%        PSF = psfGauss(dim, sigma);
%
%  This function constructs a gaussian blur PSF. 
%
%  Input: 
%    dim  -  desired dimension of the pointspread function
%            e.g., PSF = psfGauss([60,60]) creates a 60-by-60 
%            Gaussian point spread function.
%
%  Optional input parameters:
%    sigma  -  variance of the gaussian
%              Default is sigma = 2.0.
%

if ( nargin < 2 )
  sigma=2.0;
end

l = length(dim);

switch l
case 1
  x = -fix(dim(1)/2):ceil(dim(1)/2)-1;
  y = 0;
  z = 0;
case 2
  x = -fix(dim(1)/2):ceil(dim(1)/2)-1;
  y = -fix(dim(2)/2):ceil(dim(2)/2)-1;
  z = 0;
case 3
  x = -fix(dim(1)/2):ceil(dim(1)/2)-1;
  y = -fix(dim(2)/2):ceil(dim(2)/2)-1;
  z = -fix(dim(3)/2):ceil(dim(3)/2)-1;
otherwise
  error('illegal PSF dimension')
end

[X,Y,Z] = meshgrid(x,y,z);
PSF = exp( -(X.^2 + Y.^2 + Z.^2) / (2*sigma^2) ) / ( sqrt((2*pi)^l)*sigma^l );

