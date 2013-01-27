% Two-D example of Heat Eqn from Asmar's PDE book

clear all, close all
 
x = [0:0.01:1]';
y = [0:0.01:1]';
[iX,iY] = meshgrid(x,y);
nt = 50; %input(' Number of time steps = ');
dt = .1;   % Time step size.
N = 5;

% Time loop
for j = 1:nt 
  u = zeros(length(x),length(y));
  % Double sum
  for m = 1:N % index on x
    for n = 1:N % index on y
        Bmn = 1600 / (pi^2*(2*m-1)*(2*n-1));
        u = u + Bmn * sin((2*m-1)*pi*iX) .* sin((2*n-1)*pi*iY) * ...
                               exp(-( (2*m-1)^2 + (2*n-1) ) * (j*dt));
    end
  end
  mesh(x,y,u)
  axis([0 1 0 1 0 100])
  drawnow
  pause(0.1)
end

