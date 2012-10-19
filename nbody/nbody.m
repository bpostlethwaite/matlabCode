% Simple demo of Newtonian n-body problem 
% (c) 1995, Howard E. Motteler
%
% reads initial set of masses nbodyin.m
% writes final set of masses nbodyout.m
%
% x, y, z	n x 1, object co-ordinates
% M,		n x 1, object masses
% D		n x n, distance between objects
% Fx, Fy, Fz	n x n, forces between objects
% Vx, Vy, Vz	n x 1, object velocitie
%
% G		gravitational constant
% dT		time increment

% defaults
clear all
close all
k  = 1;
dT = 1;
G  = 1;

k  = 200;
% dT = input('enter dT > ');
% G  = input('enter G  > ');

load('nbodyin.mat')
nbodyin = tmp;

M  = nbodyin(:,1);
x  = nbodyin(:,2);
y  = nbodyin(:,3);
z  = nbodyin(:,4);
Vx = nbodyin(:,5);
Vy = nbodyin(:,6);
Vz = nbodyin(:,7);
[n,tmp] = size(M);
fprintf(1, 'read %d bodies\n', n);

v = 2 * [min(x),max(x),min(y),max(y)];	% set initial axes

% loop on dT steps

for i = 1:k

  % distance between all pairs of objects
  
  Fx =  ones(n,1) * x' - x * ones(1,n);
  Fy =  ones(n,1) * y' - y * ones(1,n);
  Fz =  ones(n,1) * z' - z * ones(1,n);
  
  Dsq = Fx .* Fx + Fy .* Fy + Fz .* Fz + eye(n);
  D = sqrt(Dsq);
  
  % mutual forces between all pairs of objects
  
  F =  G * (M * M') ./ Dsq;
  F =  F - diag(diag(F)); 	% set 'self attraction' to 0
  Fx = (Fx ./ D) .* F ;
  Fy = (Fy ./ D) .* F ;
  Fz = (Fz ./ D) .* F ;
  
  % net force on each body
  
  Fnet_x = Fx * ones(n,1);
  Fnet_y = Fy * ones(n,1);
  Fnet_z = Fz * ones(n,1);
  
  % change in velocity:
  
  dVx = Fnet_x * dT ./ M;
  dVy = Fnet_y * dT ./ M;
  dVz = Fnet_z * dT ./ M;
  
  Vx = Vx + dVx;
  Vy = Vy + dVy;
  Vz = Vz + dVz;
  
  % change in position
  
  x = x + Vx * dT;
  y = y + Vy * dT;
  z = z + Vz * dT;
  
  plot(x,y,'o')
  %axis(v)
  pause (.1);
  
end 

%tmp = [M,x,y,z,Vx,Vy,Vz];
%save nbodyout.mat tmp -ascii
