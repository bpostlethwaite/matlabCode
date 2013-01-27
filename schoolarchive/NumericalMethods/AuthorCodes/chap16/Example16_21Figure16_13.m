% Example 16.21 -- Figure 16.13 : Fermi-Pasta-Ulam using ODE solver rk45

clear all
close all
% parameters
 k = .00025;                  % constant step size
 tf = 500;                 % final time
 tpause = 10*k;            % record values every so often
 ipp = round(tf/tpause + eps);
    
   
% given initial conditions
 omega = 100; m = 3;
 y = [1,0,0,1/omega,0,0,1,0,0,1,0,0];
 y0 = y(:);
 e0 = hfun(0,y0);
 I = zeros(4,ipp+1);
 I(1:m,1) = .5*(y0(3*m+1:4*m).^2 + omega^2*y0(m+1:2*m).^2);
 I(m+1,1) = sum(I(1:m,1));

 % integrate using matlab's ode45
 y0 = shiftdim(y); tspan = [0 500];
 [tot,yot] = ode45(@ffun, tspan, y0);
 fprintf ('ode45 error in total energy = %e \n', hfun(tot(end),yot(end,:)')-e0 )
 % plot energies
 J = zeros(m+1,length(tot));
 J(1:m,:) = .5*( yot(:,3*m+1:4*m)'.^2 + omega^2*yot(:,m+1:2*m)'.^2);
 J(m+1,:) = sum(J(1:m,:));
 clf
 plot(tot',J)
 hold on
 axis([0 500 -0.2 1.2])
 xlabel('t')
 %legend('I_1','I_2','I_3','I')
 gtext('I_1')
 gtext('I_2')
 gtext('I_3')
 gtext('I')
 


