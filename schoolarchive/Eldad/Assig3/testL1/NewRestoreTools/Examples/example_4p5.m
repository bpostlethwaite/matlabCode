%
%  Example 4.3 from 
%
%  "Iterative Methods for Image Restoration:
%     A Matlab Object Oriented Approach"
%
%   by: K.P. Lee, J.G. Nagy and L. Perrone,
%

startup

load satellite
A = psfMatrix(PSF);
As = (A + A')/2;
P = psfPrec(As, b, 'help');
x = PMR2(As, P, b, b, 15);
figure
imshow(x,[])
