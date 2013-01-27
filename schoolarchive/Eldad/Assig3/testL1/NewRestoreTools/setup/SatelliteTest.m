clc
%
%  To set paths, you should run startup.m
%  (though of course this only needs to be done once)
%
startup

%
% The data is in NewRestoreTools/TestData/
% To load the satellite data, use:
%
disp('Loading data ...')

load satellite

disp('... done loading data')

%
% "whos" will show the data:
%    PSF = PSF kernel
%    x_true = true image
%    b = blurred and noisy image
%

%
% Construct a psfMatrix object:
%
disp('Constructing psfMatrix object ...')

A = psfMatrix(PSF);

disp('... done constructing psfMatrix object.')

%
%  If all you want to do is multiply, then
%  use, for example:
%
disp('Multiplying A*x, where x is an image array ...')

bb = A*x_true;

disp('... done multiplying A*x')
%
%  If you want to multiply by transpose, use:
%
disp('Multiplying A''*x, where x is an image array ...')

bbt = A'*x_true;

disp('... done multiplying A''*x')
%
%  Note that x_true is an image array, so the result
%  in bb is an image array of the same dimension.  If
%  you want to put images into vectors, you can do that
%  too, and the result is a vector.  For example:
%
disp('Multiplying A*x and A''*x, where x is a vector ...')

bb2 = A*x_true(:);
bbt2 = A'*x_true(:);

disp('... done multiplying A*x and A''*x.')

%
%  This returns bb2 and bbt2 as vectors of length prod(size(x_true))
%

%
%  If you want to try Julianne's  HyBR codes, do this:
%     - construct a preconditioner
%     - run the code
%  Note that her code requires the "right hand side" to be a
%  vector.  Since b is an image array, we input it into her
%  HyBR code using b(:)
%
disp('Now run Julianne''s HyBR codes')
disp('   First construct preconditioner ...')

P = new_svdPrec(A, b, 'help');

disp('   ... done constructing preconditioner.')
disp('   Now run the HyBR method ...')

x = HyBR(A, b(:), P);

disp('   ... done running HyBR method.')

disp('Display results')
figure(1), clf
subplot(2,2,1), imshow(x_true,[0,max(x_true(:))]), title('True image')
subplot(2,2,2), imshow(bb,[0,max(bb(:))]), title('Result of A*x')
subplot(2,2,3), imshow(b,[0,max(b(:))]), title('Given blurred, noisy image')
subplot(2,2,4), imshow(reshape(x,size(b)), [0,max(x(:))]), title('Restored image')

