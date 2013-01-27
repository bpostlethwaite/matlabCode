% Example 13.7 -- fig 13_5: trigonometric interpolation of square wave function

 % prepare data y
 clear all
 l = 64; % input parameter
 lenij = 1;
 fprintf(' m      err_max       err_l2        rate_max        rate_l2   \n')
 rateinf = 0; ratel2 = 0;

 for ij = 1:lenij
      
 m = 2*l; x = 0 :2*pi/m: 2*pi/m*(m-1);
 y = zeros(size(x));
 I = find ((x > 2).*(x < 4));
 y(I) = 1;

 % find dft real coefficients
 [a0,a,b] = dft1e(y);
 
 % interpolate and plot
 xx = 0:.01*pi:2*pi;
 mm = length(xx); mpi = (mm+1)/2;
 %yexact = xx; yexact(mpi+1:mm) = 2*pi - xx(mpi+1:mm);
 yexact = zeros(size(xx));
 I = find ((xx > 2).*(xx < 4));
 yexact(I) = 1;
 yapprox = dft2e(xx,a0,a,b);
 
 plot(xx,yexact,xx,yapprox)
 axis([-.5 7 -.5 1.5])
 legend('exact','approx')

 err_max = max(abs(yexact-yapprox));
 err_l2 = norm(yexact-yapprox)/sqrt(m);
 if ij > 1
     rateinf = log2(err_mo) - log2(err_max);
     ratel2 = log2(err_2o) - log2(err_l2);
 end
 fprintf('%d  %e %e %e %e \n',m,err_max,err_l2,rateinf,ratel2)
 err_mo = err_max;
 err_2o = err_l2;
 l = l*2;
 
 end