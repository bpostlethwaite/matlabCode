% Examp16.6 -- Figures 16.5, 16.6 : simple predator-prey model

y0 = [80,30]';    % initial data
tspan = [0,100];  % integration interval
h = .01;          % constant step size

[tout,yout] = rk4(@func,tspan,y0,h);

    figure(1)
    plot(tout,yout)
    xlabel('t')
    ylabel('y')
    legend('y_1','y_2')
    
    figure(2)
    plot(yout(1,:),yout(2,:))
    xlabel('y_1')
    ylabel('y_2')
   

 