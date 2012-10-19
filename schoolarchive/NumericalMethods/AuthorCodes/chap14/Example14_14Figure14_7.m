% Example 14.14 -- Figure 14.7 : solving a nonlinear partial differential
% equation with a spectral method
%
% The nonlinear Schroedinger 
%        i u_t + u_xx + gam|u|^2u = 0 
% with periodic BC is solved by splitting using spectral for u_t=iu_xx
% and exact for the nonlinear part.
% Use Strang splitting to reduce splitting error to 2nd order.

clear all
close all
h = input ('enter h = step size in x : ');
dt = input ('enter dt = time step : ');
t0 = input ('enter initial time t0 : ');
tf = input ('enter time target tf : ');
iploft = input('enter how often refresh plot (0=10) : ');
if iploft == 0, iploft = 10; end
watert = 200; waterx = 100;
nsteps = (tf-t0) / dt;

% initial conditions on the spatial mesh
  gam = 1;
  c1 = 1; c2 = 0.1;
  delta = 25;
  alfa = 0.5;
  al = -20; ar = 80;
  scale = 2*pi/(ar-al); 
  x = (al+h:h:ar)';
  n = length(x);
  if t0 == 0
    u = sqrt(2*alfa)*(exp(i*.5*c1*x).*sech(x*sqrt(alfa)) +...
         exp(i*.5*c2*(x-delta)).*sech((x-delta)*sqrt(alfa)));
  end
  

% Prepare constant factor array for advancing half a step 
expfac = exp(-i*dt/2* scale^2* [0:n/2 -n/2+1:-1]'.^2);

% Initial hamiltonian and solution norm
uc = conj(u);
e0 = sum((u(2:n)-u(1:n-1)).*(uc(2:n)-uc(1:n-1)))/(h^2) - ...
     .5*gam*sum(u(1:n).^2 .* uc(1:n).^2);
% the above line uses finite differences. Use spectral differentiation:
onedif = i*scale * [0:n/2-1 0 -n/2+1:-1]';
uh = fft(u);  wh = onedif .* uh; w = ifft(wh); wc = conj(w);
e0s = sum(w.*wc) - .5*gam*sum(u.^2 .* uc.^2);

normu0 = norm(u);
fprintf ('initial discrete Ham = %e %e\n',e0,e0s)
fprintf ('initial solution norm = %e\n',normu0)

% main loop of splitting method
for ii = 1:nsteps,

  % advance dt/2 for iu_xx
  v_hat = fft(u);
  w_hat = expfac .* v_hat;
  u = ifft(w_hat);

  % advance nonlinear term exacttly
  u = u.* exp(i*dt*gam * abs(u).^2);

  % advance dt/2 for iu_xx
  v_hat = fft(u);
  w_hat = expfac .* v_hat;
  u = ifft(w_hat);

  % plot and display discrepancy in discrete H
  if mod(ii,iploft)==0, plot(x(1:n),abs(u(1:n)));drawnow; pause(0); end
  if mod(ii,watert)==0
    % save for one display at end
    jj = ii/watert; nn = n/waterx;
    if (jj < 2000)
      uss (jj,1:nn) = abs(u(1:waterx:n)); tss(jj) = ii*dt;
      xss = x(1:waterx:n);  % save for plot
    end
  end
    
end
    
% error in hamiltonian and solution norm
uc = conj(u);
de = sum((u(2:n)-u(1:n-1)).*(uc(2:n)-uc(1:n-1)))/(h^2) - ...
     .5*gam*sum(u(1:n).^2 .* uc(1:n).^2) - e0;
% calculate err-ham also using spectral differentiation 
uh = fft(u);  wh = onedif .* uh; w = ifft(wh); wc = conj(w);
des = sum(w.*wc) - .5*gam*sum(u.^2 .* uc.^2) - e0s;

normu = norm(u) - normu0; 
fprintf ('after %d time steps, rel Ham err = %e %e\n',nsteps,de/e0,des/e0s)   
fprintf ('after %d time steps, rel norm err = %e\n',nsteps,normu/normu0)

figure(11)
clf
waterfall(xss,tss,uss)
xlabel('x')
ylabel('t')
zlabel('|\psi|')
view(-10,65)
figure(1)