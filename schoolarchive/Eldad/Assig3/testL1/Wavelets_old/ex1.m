 
    h = daubcqf(6);  load lena; 
   noisyLena = lena + 25 * randn(size(lena));
    figure; colormap(gray); imagesc(lena); title('Original Image');
     figure; colormap(gray); imagesc(noisyLena); title('Noisy Image'); 
    Denoise lena with the default method based on the DWT
  [denoisedLena,xn,opt1] = denoise(noisyLena,h);
   figure; colormap(gray); imagesc(denoisedLena); title('denoised Image');