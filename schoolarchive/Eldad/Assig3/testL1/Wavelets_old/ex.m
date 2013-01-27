x1=wavread('late_cafe');
x2=x1(1:65536);
s1=auread('late_cafe_out');
s2=s1(1:65536);
%fs=128;
%t=(1:128)/128;
%c1=sin(2*pi*t*8);
%c2=sin(2*pi*t*10);
%c3=sin(2*pi*t*60);
%c4=c1+c2+c3;
%plot(c4)
   h = daubcqf(6); 
   %[s,N] = makesig('HiSine',65536);
   %n = randn(1,N);
    %x = x2 + n/10;     % (approximately 10dB SNR)
     % figure;plot(x2);hold on;plot(s2,'r');

     %Denoise x with the default method based on the DWT
      [xd,xn,opt1] = denoise(x2,h);
     figure;plot(x2);hold on;plot(xd,'r');

    %Denoise x using the undecimated (LSI) wavelet transform
    [yd,yn,opt2] = denoise(x2,h,1,[]);
    figure;plot(x2);hold on;plot(yd,'r');
   fft1=fft(yd,65536);
   fft2=fft(x2,65536); 
   fft3=fft(s2,65536);
   w=(0:32767)/32768*(20000/2);
    figure;
    plot(w,abs(fft1(1:32768)),'c');
    hold on;
    figure;
    plot(w,abs(fft2(1:32768)),'r');
    hold on;
    plot(w,abs(fft3(1:32768)),'y');
    %legend('denoise','noisy','orignal');
    %wavwrite(yd,'clean')