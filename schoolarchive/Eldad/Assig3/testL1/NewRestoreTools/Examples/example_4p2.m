%
%  Example 4.2 from 
%
%  "Iterative Methods for Image Restoration:
%     A Matlab Object Oriented Approach"
%
%   by: K.P. Lee, J.G. Nagy and L. Perrone,
%

startup

load satellite
A = psfMatrix(PSF);
L = psfMatrix(1);
lambda = 0.00023;
Ahat = [A; lambda*L];
bhat = [b; zeros(size(b))];
P = psfPrec(Ahat, b);
x = PCGLS(Ahat, P, bhat, b, 35);
figure
imshow(x,[])
