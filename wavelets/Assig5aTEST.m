% Setup
clear all; close all;
loadtools;


%% Exercise 1: Sparse 1D Deconvolution
load source.mat
%We generate a synthetic sparse signal, with only a small number of non zero coefficients.
n = 1024;
% number of diracts
p = round(n*.03);
% location of the diracs
sel = randperm(n); sel = sel(1:p);
% signal
f = zeros(n,1);
f(sel) = 1;
f0 = f .* sign(randn(n,1)) .* (1-.3*rand(n,1));
%We load a seismic filter, which is a second derivative of a Gaussian.
% dimension of the signal
n = 1024;
% width of the filter
s = 5;
% second derivative of Gaussian
t = (-n/2:n/2-1)';
h = (1-t.^2/s^2).*exp( -(t.^2)/(2*s^2) );
h = h-mean(h);
% recenter the filter for fft use
h1 = fftshift(h);
%To ease implementation, we define a filtering function.
filter = @(u)real(ifft(fft(h1).*fft(u)));
% Display the filter and its Fourier transform.
% Fourier transform (normalized)
hf = real(fftshift(fft(h1))) / sqrt(n);
% display
q = 200;
clf;
subplot(2,1,1);
plot(-n/2+1:n/2, h);
axis([-q q min(h) max(h)]);
title('Filter, Spacial (zoom)');
subplot(2,1,2);
plot(-n/2+1:n/2, hf);
axis([-q q 0 max(hf)]);
title('Filter, Fourier (zoom)');

%We compute blurry noisy measurements y=h*f0+w where w is the noise,.
% noise level
sigma = .06;
% noise
w = sigma*randn(n,1);
% blurring + noise
y = filter(f0) + w;
%Display signals and measurements.
clf;
subplot(2,1,1);
plot_sparse_diracs(f0); axis('tight');
title('Signal f0');
subplot(2,1,2);
plot(y); axis('tight');
title('Measurements y=h*f0+w');

% gradient descent step size
tau = 1.95 / max(abs(fft(h)))^2;
% regularization strenght
lambda = .1;
% initialization
fSp = y;

%For such an ill-posed inverse problem, the number of iterations needs to be very large.
niter = 2000;
%The first step, like for pseudo inverse computations is a gradient descent step.
clf
h = plot(fSp,'YDataSource','fSp');
for ii = 1:niter
    fSp = fSp - tau*filter( filter(fSp) - y );
    %The second step is a soft thresholding. The threshold must be set to lambda*tau.
    fSp = perform_thresholding(fSp, lambda*tau, 'soft');
    %refreshdata(fSp,'caller') % Evaluate y in the function workspace
	%drawnow;
end
