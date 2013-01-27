% Example 16.20 -- Figure 6.11

clear all
close all
y0 = [.994,0,0,-2.00158510637908252240537862224]';  % initial data
isw = 1;

if isw == 2  % use ode45
tspan = [0 17.1];
[tout,yout] = ode45(@func,tspan,y0);
yout = yout';
figure(isw)
plot(yout(1,:),yout(3,:))
    xlabel('u_1')
    ylabel('u_2')

else % use RK4 with sppecified h
    
tspan = [0,17.1];
h = tspan(2)/10000;
[tout,yout] = rk4(@func,tspan,y0,h);

    figure(1)
    plot(tout,yout(1,:),tout,yout(3,:))
    xlabel('t')
    legend('u_1','u_2')
    
    figure(2)
    plot(yout(1,:),yout(3,:))
    xlabel('u_1')
    ylabel('u_2')
   
end
 