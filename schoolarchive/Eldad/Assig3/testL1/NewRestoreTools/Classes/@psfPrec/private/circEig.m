function E = circEig(PSF, center)
%
%          E = circEig(PSF, center);
%
%  Compute the eigenvalues of the "Strang" circulant preconditioner,
%  that is, find the eigenvalues of the circulant matrix that minimizes
%
%                || A - C ||
%                           1
%
%  where A is a psfMatrix.
%
%  Input: 
%      PSF  -  array containing a single PSF
%   center  -  array containing the indices of the center of the PSF
%
%  Output:
%        E  -  array containing the complex eigenvalues of the
%              circulant preconditioner
%

%  J. Nagy  1/2/02

PSF = PSF - min( PSF(:) );
%PSF = PSF / sum( PSF(:) );
E = fftn( circshift(PSF, -(center - 1)) );
