%  Figure 11.8 : plot B-splines
%
clear all
close all
   k = 3; % cubic
   m = 2; % C^m
   nt = 6;
   t = 0:1:nt;  % abscissae
   
   % knots
   x(1) = t(1);
   for j = 2:k, x(j) = t(1)-j*1.e-5; end
   for i = 2:nt
     for j = 1:k-m
       x = [x,t(i)];
     end
   end  
   for j = 1:k, x = [x,t(end)+j*1.e-5]; end  
   
   xx = 0:.01:nt;  % evaluation points for plotting
   start = 1; fin = 8;
   for j=start:fin % the jth B-spline
     tt = x(j:j+k);
     for l=1:length(xx)
       for ll = 1:length(tt)
         yy(ll) = max(xx(l)-tt(ll),0)^(k-1);
       end
       [coef,table] = divdif(tt,yy);
       Bj(l) = (-1)^(k)*(x(j+k)-x(j))*table(k+1,k+1);
     end

     plot(xx,Bj)
     if j==start, hold, xlabel('x'), ylabel ('B_j'), end
     if j == ceil((start+fin)/2)
       plot(xx,Bj,'LineWidth',2,'Color',[.6,0,0])
     end

   end
