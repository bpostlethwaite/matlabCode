%
%  This shows how to use the iterative method HyBR with the
%  RestoreTools package.   Currently to use HyBR, you need
%  to vectorize the blurred image, and the result you compute
%  will be a vector as well.  That is:
%     * Use b(:) instead of b when calling HyBR
%     * the solution, x, will be a vector.  To display,
%       use reshape(x, m, n), where size(image) = m-by-n
%  Some examples below illustrate how to use HyBR.  For more 
%  information, see help HyBR.
%
%  Reference: "A modified GCV method for Lanczos hybrid regularization"
%              J. Chung, J. Nagy and D. O'Leary
%

%
% Easy example is the satellite:
%    * Use zero boundary conditions
%    * Use new_svdPrec as a preconditioner, and allow "help"
%      to choose appropriate regularization for it.
%
clc
disp('Satellite example')
load satellite
A = psfMatrix(PSF);
P = new_svdPrec(A, b, 'help');
[x, exit] = HyBR(A, b(:), P);
figure(1), imshow(b,[]), title('blurred satellite')
figure(2), imshow(reshape(x, size(b)), [0,max(x(:))]), title('reconstructed satellite')
disp(sprintf('number of HyBR iterations = %d', exit.iterations))

%
% In the case of the grain data, reflexive boundary conditions
% are most appropriate:
%
disp(' ')
disp('Grain example')
load Grain
A = psfMatrix(PSF, 'reflexive');
P = new_svdPrec(A, b, 'help');
[x, exit] = HyBR(A, b(:), P);
figure(3), imshow(b,[]), title('blurred grain')
figure(4), imshow(reshape(x, size(b)), [0,max(x(:))]), title('reconstructed grain')
disp(sprintf('number of HyBR iterations = %d', exit.iterations))

%
% In the case the star_cluster data, we need to do two things:
%   * First, HyBR has trouble determining the size of the "blurring matrix"
%     from the psfMatrix if the PSF is not the same size as the blurred
%     image.  To get around this, we use padarray to first increase the
%     size of the PSF.
%   * Second, this problem is not ill-conditioned, so asking for "help"
%     when constructing new_svdPrec will give a preconditioner that is not
%     optimal.  So here we set the regularization in the preconditioning to
%     0.
%   * Note also that we rescale the data for display.
%
disp(' ')
disp('star cluster example')
load star_cluster
PSF = padarray(PSF1, size(b)-size(PSF1), 'post');
A = psfMatrix(PSF);
P = new_svdPrec(A, b, 0);
[x, exit] = HyBR(A, b(:), P);
figure(5), imshow(b,[50,500]), title('blurred star cluster')
figure(6), imshow(reshape(x, size(b)), [50,500]), title('reconstructed star cluster')
disp(sprintf('number of HyBR iterations = %d', exit.iterations))