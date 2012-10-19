% Example 13.3 -- Figure 13.3 : trigonometric interpolation of hat function

clear all
close all
% prepare data y
 l = 2; % input parameter
 m = 2*l;
 for ii = 1:m
   x(ii) = 2*pi*(ii-1)/m;
   y(ii) = x(ii);
   if x(ii) > pi
      y(ii) = 2*pi - x(ii);
   end
 end

 % find dft real coefficients
 [a0,a,b] = dft1e(y)
 
 % interpolate and plot
 xx = 0:.01*pi:2*pi;
 mm = length(xx); mpi = (mm+1)/2;
 yexact = xx; yexact(mpi+1:mm) = 2*pi - xx(mpi+1:mm);
 yapprox = dft2e(xx,a0,a,b);
 
 plot(xx,yexact,xx,yapprox)
 axis([-.5 7 -.5 3.5])
 xlabel('x')
 legend('exact','approx')

 % err_max = max(abs(yexact-yapprox))