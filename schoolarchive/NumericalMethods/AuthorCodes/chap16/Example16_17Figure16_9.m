% Example 16.17 -- Figure 16.9 : Program to solve a simple heat eqn in 
% two space variables on the unit square.

clear all
close all
tf = .1;    % final time
h = .01;   % space step size
yy = [0:h:1];  % for plotting
J = length(yy) - 2;
Jp1 = J + 1; Jp2 = J+2;
u0 = 25 * ones(Jp2,Jp2);   % initial conditions
Lap = -delsq(numgrid('S',Jp2));
isw = 1;

if isw == 1                % Forward Euler in time
  k = 1/12*h*h;
  nf = tf / k;
  mu = k /h^2;
  A = speye(J^2) + mu * Lap;
  v = reshape(u0(2:Jp1,2:Jp1),J^2,1);
  
  for n=1:nf
    v = A*v;
  end

elseif isw == 2         % Crank-Nicolson, wastefully
  k = h;
  nf = tf / k;
  mu = k /h^2;   
  Al = speye(J^2) - .5*mu*Lap; Ar = speye(J^2) + .5*mu*Lap;
  v = reshape(u0(2:Jp1,2:Jp1),J^2,1);
  
  for n=1:nf
    v = Al \ (Ar*v);
  end
  
end

vn = u0 * 0;
vn(2:Jp1,2:Jp1) = reshape(v,J,J);
%plot solution
xx = [0:h:1]';
mesh(xx,yy,vn)
xlabel('y')
ylabel('x')
zlabel('u')
