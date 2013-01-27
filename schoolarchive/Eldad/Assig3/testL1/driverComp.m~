%driver for the inversion
clear all
close all

%%
disp(' ===== Preparing inverse problem =======');
fprintf('\n\n');

%%%%%%%
% Jim's package for the deblurring problem
addpath NewRestoreTools/setup
startupNRT;

% Wright's package for wavelets and L1
addpath GPSR_6.0

% wavelets package for L1 basis
addpath Wavelets

% TV
addpath TV


%%
%%%  Get the matrix and data %%%%%%%%%%%%%
example = 2;
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
    PSF = getPSFe(256,1e-2);
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
subplot(2,2,2)
imagesc(reshape(b,size(PSF,1),size(PSF,2))); colormap(gray);
title('The data');
axis off square tight
drawnow;

subplot(2,2,1);
imagesc(x_true); colormap(gray);
title('True image');
axis off square tight
drawnow;

%%
%%%%% setup the wavelet basis
%
wav = daubcqf(2);
W = @(x) midwt(x,wav,3);
WT = @(x) mdwt(x,wav,3);

%%
%%%%%% Now apply different solvers %%%%%

%% Jim's solver

fprintf('========================\n\n');
disp('Using L2 type solver')

fprintf('========================\n\n\n\n');

tt1 = cputime;
P = new_svdPrec(A, reshape(b,n1,n2), 'help');
xHyBR = HyBR(A, b(:));

xHyBR(find(xHyBR<0)) = 0;
tt1 =  cputime - tt1;

subplot(2,2,3);
imagesc(reshape(xHyBR,n1,n2))
title('HyBR');
colormap(gray)
axis off square tight
drawnow;

%% Steve's solver
fprintf('========================\n\n');
disp('GPSR_6.0 ');
fprintf('========================\n\n\n\n');

% move from a vector to image
rshp1 = @(x) (x(:));
rshp2 = @(x) (reshape(x,n1,n2));

AD = @(x) rshp2(A*( rshp1(W(rshp2(x) )  )));
ADT = @(x) rshp2(WT(rshp2(A'*rshp1(x) ) ));
%
% regularization parameter (image dependant)
regpar = 1e-6; % for satelite
%regpar = 1e0; % for MRI
%regpar = 1e0; % for Barbara

% set tolA
tolA = 1.e-3;
%
tt2 = cputime;
%
b2 = rshp2(b);
% [theta,theta_debias] = GPSR_Basic(b2, AD, tau,...
% 	'Debias',0,...
% 	'AT',ADT,... 
%     'True_x',WT(x_true),...
% 	'Initialization',ADT(b2),...
% 	'StopCriterion',4);

theta = IST(b2,AD,regpar,...
                	'Debias',0,...
	                'AT',ADT,... 
                   'True_x',WT(x_true),...
	               'Initialization',ADT(b2),...
	               'StopCriterion',1,...
	               'ToleranceA',tolA);


tt2 = cputime - tt2;
%
xGPSR = reshape(W(theta),n1,n2);

subplot(2,2,4)
imagesc(xGPSR)
title('GPSR_{6.0}');
colormap(gray)
axis off square tight
drawnow;


%% TV
tt3 = cputime;
maxiter = 150;
%regpar  = 1e-3;  % for satelite
regpar  = 3e-2;  % for MRI
%regpar  = 1e0;  % for  Barbara

xTV = TVsolve(A,b,regpar,[n1,n2],maxiter);

subplot(2,2,2);
imagesc(reshape(xTV,n1,n2));
title('TV');
colormap(gray)
axis off square tight
drawnow;
tt3 = cputime - tt3;



%%
%%% Some numbers for comparison
t1 = norm(x_true(:) - xHyBR(:))/norm(x_true(:));
t2 = norm(x_true(:) - xGPSR(:))/norm(x_true(:));
t3 = norm(x_true(:) - xTV(:))/norm(x_true(:));

mis1 = norm(A*xHyBR(:) - b(:))/norm(b(:));
mis2 = norm(A*xGPSR(:) - b(:))/norm(b(:));
mis3 = norm(A*xTV(:) - b(:))/norm(b(:));


fprintf('                   HyBR         GPSR_6.0          TV\n');
fprintf('misfit              %3.2e     %3.2e     %3.2e\n',mis1,mis2,mis3);
fprintf('||x-xtrue||       %3.2e     %3.2e     %3.2e\n',t1,t2,t3);
fprintf('time                %3.2e     %3.2e     %3.2e\n',tt1,tt2,tt3);


return