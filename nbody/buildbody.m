% Generate a set of "bodies" for the n-body program
%
% output is an n by 7 ascii matrix nbodyin.mat
%
% bodies are represented as an n x 7 ascii matrix with
% columns: mass, x, y, z position, x, y, z velocity


m = 1; % mass
n  = 3; % number of bodies

% random initial mass, position, and velocity

M  =  m * (rand(n,1) * 2 + .5);

x  = 20 * (rand(n,1)-0.5);
y  = 20 * (rand(n,1)-0.5);
z  = zeros(n,1);

Vx = (rand(n,1)-0.5) * 2;
Vy = (rand(n,1)-0.5) * 2;
Vz = zeros(n,1);

tmp = [M,x,y,z,Vx,Vy,Vz];
save('nbodyin.mat','tmp')