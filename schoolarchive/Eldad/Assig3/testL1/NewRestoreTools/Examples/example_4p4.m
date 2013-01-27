%
%  Example 4.4 from 
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
x = PMRNSD(A, P, b, b, 5);
figure
imshow(x,[])
