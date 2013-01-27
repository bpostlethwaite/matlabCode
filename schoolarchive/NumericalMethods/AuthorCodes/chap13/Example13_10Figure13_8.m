% Example 13.10 -- Figure 13.8 : trigonometric interpolation of log(1+x)

  % prepare data y
  clear all
  close all
  l = 16; % input parameter
  m = 2*l; x = .5*pi/m :pi/m: (m-.5)*pi/m;
  t = 2*x;
  y = log(t+1);
 
  % find dft real coefficients
  [a] = dct1e(y,2);
 
  % interpolate and plot
  xx = .5*pi/m:.001*pi:pi-.5*pi/m; tt = 2*xx;
  yexact = log(tt+1);
  yapprox = dct2e(xx,a,2);
 
  plot(tt,yexact,'--',tt,yapprox,'LineWidth',2)
  axis([-.5 7 -.5 2.5])
  xlabel('x')
  ylabel('log (1+x)')
  legend('exact','approx')

  err_max = max(abs(yexact-yapprox))
  
