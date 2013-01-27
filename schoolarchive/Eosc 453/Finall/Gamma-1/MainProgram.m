%% --------------------------- Main Program -------------------------------
% Assignment 3, EOS 453 -- Eric Deal, Charlie Beard, Saman Monfared
close all
clear all
for ii = 1:5
global Rbo awo gamma qa dt Rknow  t_length 

% set parameters
T1 = 270; % temperature at bottom of crust
Tn = 103; % temperature at surface
Rbo = 1560000; % radius of ocean
Rc = 1565000; % radius of planet
n = (Rc - Rbo)/100; % number of intervals in ice crust
To = linspace(T1,Tn,n); % initial conditions
awo = (ii)*1e-7; % heat production per unit vol water due to tidal dissapation
qa = 0.01; % heat flux from core to ocean
gamma = -1;



% calculate time step
t_length = 1e6; % length of simulation in years
yr_s = 31557600; % number of seconds in a year
timespan = linspace(0,t_length*yr_s); % Calculate time span
% choose a time step for the ice formation based on ODE solver in seconds
year = 31557600;
month = 2629800;
week = 606880;
day = 86400;
timestep = 100*year;
dt = (t_length*yr_s*(1/timestep) + 50)\(t_length*yr_s);
tic
options  = odeset('maxstep',timestep);
[t,T] = ode45(@ode,timespan,To,options); % pass ODE to ODE solver
toc

t_yr = t/yr_s; % number of years in time span


%% ---------------------------- Plot Results ------------------------------

clim=[100 273];
clim2 = [100 273];
depth = linspace(0,(Rc-Rbo),n)/1000;
time = linspace(0,t_length);
year = 6000;

%temperature limitations of plot
figure(ii);
imagesc(time,depth,flipud(real(T)'),clim);
colorbar;
h=colorbar ();
xlabel(h,'Temp(K)');
xlabel ('time (years)');
ylabel ('depth of ice crust (km)');
title (sprintf('temperature evolution of ice crust of europa-- gamma = %3.1f -- aw = %2.2f',gamma,awo));

Rtime = linspace(0,t_length,length(Rknow));

figure(10+ii)
plot(Rtime,Rc - Rknow,'k')
xlabel('time (yrs)')
ylabel(sprintf('Thickness of crust (m) -- aw = %2.2f',awo))


[Rt] = tpole(T);
TT{ii} = T;
tempprofile(ii) = Rt;
v = size(T); 
figure(20+ii)
plot(T(1,:),fliplr(depth),'r',T(v(1),:),fliplr(depth),'g.')
set(gca,'ydir','reverse')
xlim ([100 350])
xlabel ('temp (k)')
ylabel ('depth (km)')
title (sprintf('temperature-depth snapshots of europa ice crust -- gamma = %3.1f -- Rt = %2.3f -- aw = %2.2f',gamma,Rt,awo));
legend ('time = 0 yr', sprintf('time = %i Kyr',(t_length)/1000))
filename = sprintf('gamma_%3.0f_aw_%2.0f',gamma,awo);
save(filename,'T','Rt');
end
tempprofile