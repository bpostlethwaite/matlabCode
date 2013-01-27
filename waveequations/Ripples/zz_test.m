clear all; close all;
c=1; % speed of wave;
dx=1; % space step;
dt=0.1; % time step;
szx=200;
szy=200; % sizes
tm=1000; % time
k=0.002; % decay factor

<<<<<<< HEAD
dp=0.01; % droplet probability per one time sterp
dsz=3; % droplet size
da=20; % droplet amplitude

=======
>>>>>>> origin/master
x=0:dx:szx;
y=0:dx:szy; % space
t=0:dt:tm; % time

Lx=length(x);
Ly=length(y);

u=zeros(Ly,Lx); % initial value
<<<<<<< HEAD
=======
u(50,70)=100; 
>>>>>>> origin/master
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
<<<<<<< HEAD


% droplet as gaussian 
xd=-2*dsz:dx:2*dsz;
yd=-2*dsz:dx:2*dsz;
[Xd,Yd] = meshgrid(xd,yd);
Zd=-da*exp(-(Xd/dsz).^2-(Yd/dsz).^2);
point = [0,0];   % This is the current clicked cursor position in figure
point0 = point;  % This is the previous clicked cursor position

while true
=======
for tt=t
>>>>>>> origin/master
    un=(2-kdt)*u+(kdt-1)*uo+c1*conv2(u,D,'same'); % new
    
    uo=u; % curent become old
    u=un; % new become current
    
    if mod(lc-1,dlc)==0
        set(hi,'CData',u);
        drawnow;
    end
    
<<<<<<< HEAD
    
    if mod(lc-1,1)==0  % Don't want to check too much for speed
        pointg = get(ha,'CurrentPoint'); % Get current clicked point
        point = [pointg(1,1),pointg(1,2)];
        if point ~= point0 % Compare with previous, if not equal, add droplet
            if point(1) > 0 && point(1) < Lx && point(2) > 0 && point(2) < Ly
                %[point(1),point(2)]
                point0 = point;
                x0d=2*dsz+1+(szx-4*dsz-2)*point(1)/Lx;
                y0d=2*dsz+1+(szy-4*dsz-2)*point(2)/Ly; % droplet center
                x0d=round(x0d);
                y0d=round(y0d);
                u(y0d-2*dsz:y0d+2*dsz,x0d-2*dsz:x0d+2*dsz)=...
                    u(y0d-2*dsz:y0d+2*dsz,x0d-2*dsz:x0d+2*dsz)+Zd;
            end
        end
    end
    
    
=======
>>>>>>> origin/master
    lc=lc+1;
end


