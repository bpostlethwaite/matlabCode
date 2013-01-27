%
%  Example 4.1 from 
%
%  "Iterative Methods for Image Restoration:
%     A Matlab Object Oriented Approach"
%
%   by: K.P. Lee, J.G. Nagy and L. Perrone,
%

startup

load satellite
A = psfMatrix(PSF);
P = psfPrec(A, b, 'help');
x = PCGLS(A, P, b, b, 15);
figure
imshow(x,[])
