% Example 13.5 : extending Example 13.4 with r=10 and 2. record errors

 clear all
 r = 10; 
 l = 16; lenij = 3; % input parameters
 
 fprintf(' m      err_max       err_l2        rate_max        rate_l2   \n')
 rateinf = 0; ratel2 = 0;

 for ij = 1:lenij
     
   m = 2*l; pi2 = 2*pi;
   x = 0 : pi2/m : (m-1)*pi2/m; % abscissae on [0,pi2]
   t = 3/pi2*x - 1;             % abscissae on [-1,2]
   y = t.^2 .* (t+1).^r .* (t-2).^r - ...
     exp(-t.^2) .* sin(t+1).^r .* sin(t-2).^r; % data y

   % find dft real coefficients
   [a0,a,b] = dft1e(y);
 
   % interpolate on fine mesh, plot and find max error
   xx = 0:.01*pi:2*pi;
   tt = 3/pi2*xx - 1;
   yexact = tt.^2 .* (tt+1).^r .* (tt-2).^r - ...
     exp(-tt.^2) .* sin(tt+1).^r .* sin(tt-2).^r;
   yapprox = dft2e(xx,a0,a,b);
 
   err_max = max(abs(yexact-yapprox)./(abs(yexact)+1));
   err_l2 = norm((yexact-yapprox)./(abs(yexact)+1))/sqrt(m);
   if ij > 1
     rateinf = log2(err_mo) - log2(err_max);
     ratel2 = log2(err_2o) - log2(err_l2);
   end
   fprintf('%d  %e %e %e %e \n',m,err_max,err_l2,rateinf,ratel2)
   err_mo = err_max;
   err_2o = err_l2;
   l = l*2;
 
 end