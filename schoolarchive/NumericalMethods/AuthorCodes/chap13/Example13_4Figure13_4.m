% Example 13.4 -- Figure 13.4 : trigonometric interpolation of the function
% g(t) = t^2(t+1)^2(t-2)^2 - e^{-t^2}\sin^2(t+1)\sin^2(t-2)
%        on the interval [-1,2].

% prepare data y
 clear all
 close all
 l = 1;
 for j=1:4
   l = l*2;  
   m = 2*l; pi2 = 2*pi;
   x = 0 : pi2/m : (m-1)*pi2/m; % abscissae on [0,pi2]
   t = 3/pi2*x - 1;             % abscissae on [-1,2]
   y = t.^2 .* (t+1).^2 .* (t-2).^2 - exp(-t.^2) .* sin(t+1).^2 .* sin(t-2).^2;

   % find dft real coefficients
   [a0,a,b] = dft1e(y);
 
   % interpolate on fine mesh, plot and find max error
   xx = 0:.01*pi:2*pi; tt = 3/pi2*xx - 1;
   yexact = tt.^2 .* (tt+1).^2 .* (tt-2).^2 - exp(-tt.^2) .* sin(tt+1).^2 .* sin(tt-2).^2;
   yapprox = dft2e(xx,a0,a,b);
 
   % display results
   subplot(2,2,j)
   plot(tt,yexact,'b--',tt,yapprox)
   axis([-1 2 -2 6])
   legend('exact','approx')
   xlabel('t')

   err_max = max(abs(yexact-yapprox)./(abs(yexact)+1))
 end