% Example 8.6 -- Figures 8.5 and 8.6 : image compression with SVD

close all
colormap('gray')
load clown.mat;
figure(1)
image(X);

[U,S,V] = svd(X);
figure(2)
r = 20;
colormap('gray')
image(U(:,1:r)*S(1:r,1:r)*V(:,1:r)');