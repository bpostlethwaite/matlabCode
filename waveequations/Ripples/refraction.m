c=1; % speed of wave;
dx=1; % space step;
dt=0.1; % time step;

I=imread('matlab_logo.png');
[szy szx ch]=size(I);

%szx=200;
%szy=200; % sizes
tm=3000; % time
k=0.002; % decay factor
dp=0.01; % droplet probability per one time sterp
dsz=3; % droplet size
da=1; % droplet amplitude

gf=0.3; % gradient factor

x=0:dx:(szx-dx);
y=0:dx:(szy-dx); % space
t=0:dt:tm; % time

[X,Y] = meshgrid(x,y);

Lx=length(x);
Ly=length(y);

u=zeros(Ly,Lx); % initial value
%u(50,70)=100; 
uo=u; % previose = curent => velocties =0

close all;
hf=figure;
ha=axes;
hi=imshow(I);
% hi=imagesc(x,y,u);
% colorbar;
% set(ha,'clim',[-1 1]);
% axis equal;
% 
% % cyan colormap:
% bl1=linspace(0,1,64)';
% bl2=linspace(0.0,1,64)';
% bl3=linspace(0.0,0,64)';
% bl=[bl3 bl2 bl1];
% colormap(bl);

D=[0 1 0; 1 -4 1; 0 1 0]; % 2d laplace operator

% sobel operator:
Gx=[-1 0 1; 
    -2 0 2;
    -1 0 1];
Gy=[-1 -2 -1;
    0 0 0;
    1 2 1];

kdt=k*dt;
c1=dt^2*c^2/dx^2;
lc=1;
dlc=15;
Ii=I;

% droplet as gaussian 
xd=-2*dsz:dx:2*dsz;
yd=-2*dsz:dx:2*dsz;
[Xd,Yd] = meshgrid(xd,yd);
Zd=-da*exp(-(Xd/dsz).^2-(Yd/dsz).^2);

for tt=t
    un=(2-kdt)*u+(kdt-1)*uo+c1*conv2(u,D,'same'); % new
    
    uo=u; % curent become old
    u=un; % new become current
    
    if mod(lc-1,dlc)==0
        Xi=X+gf*conv2(u,Gx,'same');
        Yi=Y+gf*conv2(u,Gy,'same');
        Ii(:,:,1) = interp2(X,Y,double(I(:,:,1)),Xi,Yi,'*linear',255);
        Ii(:,:,2) = interp2(X,Y,double(I(:,:,2)),Xi,Yi,'*linear',255);
        Ii(:,:,3) = interp2(X,Y,double(I(:,:,3)),Xi,Yi,'*linear',255);
        set(hi,'CData',Ii);
        %set(hi,'Cdata',u);
        drawnow;
    end
    
    % droplets:
    if rand<dp
        x0d=2*dsz+1+(szx-4*dsz-2)*rand;
        y0d=2*dsz+1+(szy-4*dsz-2)*rand; % droplet center
        x0d=round(x0d);
        y0d=round(y0d);
        u(y0d-2*dsz:y0d+2*dsz,x0d-2*dsz:x0d+2*dsz)=...
            u(y0d-2*dsz:y0d+2*dsz,x0d-2*dsz:x0d+2*dsz)+Zd;
    end
    
    
    
    lc=lc+1;
end


