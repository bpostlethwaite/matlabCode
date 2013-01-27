tmin=0; tmax=100;
x=linspace(tmin,tmax,100);
mag = 3;
a=2;
ntime=100;
itime=x;
mod=ones(1,ntime);
n=4;
for ii = 1:n
    
mod(1,ii*ntime/n-round(0.5*ntime/n)-0.5*a:ii*ntime/n-round(0.5*ntime/n))=3;

end



