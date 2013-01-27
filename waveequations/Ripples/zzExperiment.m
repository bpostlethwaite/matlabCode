clear all; close all;
<<<<<<< HEAD
c=1; % speed of wave;
dx=1; % space step;
dt=0.1; % time step;
=======
c = 1; % speed of wave;
dx = 1; % space step;
dt = 0.1; % time step;
>>>>>>> origin/master
szx=200;
szy=200; % sizes
tm=1000; % time
k=0.002; % decay factor

x=0:dx:szx;
y=0:dx:szy; % space
%t=0:dt:tm; % time % Going to make infinite loop

Lx=length(x);
Ly=length(y);

u=zeros(Ly,Lx); % initial value
u(50,70)=100;
uo=u; % previose = curent => velocties =0

close all;
hf=figure;
ha=axes;
hi=imagesc(x,y,u);
colorbar;
set(ha,'clim',[-1 1]);
axis equal;

% cyan colormap:
bl1=linspace(0,1,64)';
bl2=linspace(0.0,1,64)';
bl3=linspace(0.0,0,64)';
bl=[bl3 bl2 bl1];
colormap(bl);

D=[0 1 0; 1 -4 1; 0 1 0]; % 2d laplace operator

kdt=k*dt;
c1=dt^2*c^2/dx^2;
lc=1;
dlc=15;


while 1
    un=(2-kdt)*u+(kdt-1)*uo+c1*conv2(u,D,'same'); % new
 
    uo=u; % curent become old
    u=un; % new become current
    
    if mod(lc-1,dlc)==0
        set(hi,'CData',u);
        drawnow;
    end
    
    lc=lc+1;
end


