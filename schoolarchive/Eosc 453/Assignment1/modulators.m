clear all


yrmin=1800; yrmax=2150; npoints=500;
yr=linspace(yrmin,yrmax,npoints);
interval=yrmax-yrmin;
time=linspace(0,interval,npoints);

y=cos(time/12)/4+.8;

plot(yr,y);



