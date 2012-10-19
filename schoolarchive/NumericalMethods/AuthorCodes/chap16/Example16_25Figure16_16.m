% Example 16.25 -- Figure 16.16 : waterfall plot for wave equation
% u_tt = u_xx 
%         -10 < x < 10,  u(0,x) = exp(-alfa x^2), u_t(0,x) = 0
% Dirichlet BC

clear all
figure(1)

% default values
alfa = 1;
J = 500; % number of spatial mesh points
dt = .02; % time step
% or read them in
alfa = input(' enter alfa for initial data : ');
J = input ('enter number of spatial mesh points J : ');
dt = input ('enter time step : ');


% Grid and initial data:
a = -10; b = 10;
t = 0;  tf = 40; N = round(tf/dt);
h = (b-a)/J; x = a:h:b;
  
v0 = exp(-alfa*x.^2); v0(1) = 0; v0(J+1) = 0;
vm1 = v0;

% Time-stepping by FD leap frog formula:
     
  mu = dt/h;
  vnew(2:J) = v0(2:J) + .5*mu^2* ...
          (v0(3:J+1)-2*v0(2:J)+v0(1:J-1)); % time step -1 
  vnew(1) = 0; vnew(J+1) = 0; 
  vm1 = vnew; 
  v = v0; vold = vm1;
  
  for n = 1:N 
      
      vnew(2:J) = 2*v(2:J) - vold(2:J) + mu^2* ...
          (v(3:J+1)-2*v(2:J)+v(1:J-1));
      vold = v; v = vnew;
      
      if mod(n,10)==0,
        plot(x(2:J+1),vnew(2:J+1)); 
        axis([-10,10,-1.1,1.1]); drawnow; pause(0);
        % save for one display at end
        ii = n/10; nn = J/10;
        if (ii < 1000)
          uss (ii,1:nn) = vnew(2:10:J+1); tss(ii) = n*dt;
          xss = x(2:10:J+1);  % save for plot
        end
      end

  end
  %plot(x,vnew(2:J+1),'b')
  %xlabel x, ylabel u
  %err = norm(vnew-v0,inf)
  
figure(3)
clf
waterfall(xss,tss,uss)
xlabel('x')
ylabel('t')
zlabel('u')
view(-10,65)
  
