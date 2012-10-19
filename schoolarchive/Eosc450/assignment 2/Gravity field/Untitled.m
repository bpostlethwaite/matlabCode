load grav67.doc
x=1;
while x<=1173
f(x,1)=grav67(x,2);
l(x,1)=grav67(x,3);
DgF(x,1)=grav67(x,5);
f1(x,1)=f(x,1)*pi/180;
x=x+1;
end
for x=1:1173
    go67(x,1)=(9.78031846*(1+(0.0053024*sin(f1(x,1))^2)-(0.0000059*sin(2*f1(x,1))^2)))*(10^5);
    go80(x,1)=(9.78032677*((1+(0.001931851353*sin(f1(x,1))^2))/(1-(0.0066943800229*sin(f1(x,1))^2))^(1/2)))*(10^5);
    Dg80(x,1)=DgF(x,1)+go67(x,1)-go80(x,1);
end
[x,y]=meshgrid(-120:(1/12):-114,56:(1/12):60);
z=griddata(f,l,Dg80,y,x,'v4');
load grav80.txt
for i=1:1173
    EIGEN(i,1)=Dg80(i,1)-grav80(i,5);
    
end
o=griddata(f,l,EIGEN,y,x,'v4');
contourf(x,y,z),colorbar,title('GRS80(1)')
figure
contourf(x,y,o),colorbar,title('EIGEN(1)')
load grav80rtm.txt
for i=1:1173
    MRT(i,1)=EIGEN(i,1)-grav80rtm(i,5);
end
figure
lo=griddata(f,l,MRT,y,x,'v4');
contourf(x,y,lo),colorbar


