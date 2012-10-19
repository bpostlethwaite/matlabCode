% Example 11.7 -- Figures 11.9 and 11.10

clear all
xi = [-1,0,1,0,1];
yi = [0,1,0.5,0,-1];
figure(1)
clf
plot (xi,yi,'go',xi,yi)
axis([-1.5,1.5,-1.5,1.5])
xlabel('x')
ylabel('y')
  
figure(2)
clf
ti = 0:0.25:1;
coefx = divdif(ti,xi);
coefy = divdif(ti,yi);
tau = 0:.01:1;
xx = evalnewt(tau,ti,coefx);
yy = evalnewt(tau,ti,coefy);
plot (xi,yi,'go',xx,yy)
axis([-1.5,1.5,-1.5,1.5])
xlabel('x')
ylabel('y') 