%%%%%%%%%%%%%%%%%%%%%%%%%% test TV %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Jim's package for the deblurring problem
addpath ../
addpath ../NewRestoreTools/setup
addpath ../GPSR_6.0 
startupNRT;


%%%  Get the matrix and data %%%%%%%%%%%%%
example = 4;
if example == 1,  % satellite
   load satellite
elseif example == 2 % MRI
    %load PSF128 % for the PSF
    PSF = getPSFe(128,1e-2);
    load mri
    x_true = double(D(:,:,16));
elseif example == 3,
    %load PSF128 % for the PSF
    PSF = getPSFe(128,1e-2);
    load woman2
    x_true = double(X);
elseif example == 4,
    %load PSF128 % for the PSF
    PSF = getPSFe(256,1e-3);
    X = double(imread('Camera.tif'));
    %X = 0.25*(X(1:2:end,1:2:end) + X(1:2:end,2:2:end) + ...
    %                X(2:2:end,1:2:end) + X(2:2:end,2:2:end));     
    x_true = X;
end

[n1,n2] = size(x_true);

% generate the forward modeling matrix
A = psfMatrix(PSF);
b = A*x_true(:);

% add noise
sigma = mean(abs(b(:)))/100;
noise = randn(size(b)) * sigma;
b = b+noise;

figure(1);
subplot(2,2,1)
imagesc(reshape(b,size(PSF,1),size(PSF,2))); colormap(gray);
title('The data');
axis off square tight
drawnow;

subplot(2,2,2);
imagesc(x_true); colormap(gray);
title('True image');
axis off square tight
drawnow;

%% Solve using TV
maxiter = 10;
regpar  = 0.1; 
[x] = TVsolve(A,b,regpar,[n1,n2],maxiter);