%% GravityDriver    

%  Note to Eldad: I am using my own W function that targets specific
%  depths, it is tuned to the depth of the density anomaly.
% A helpful geologist I suppose... anyway, you can change W back
% to your W  on line 53, just make it WtW = W'*W instead of the Wz's

clear all
close all

%% Define parameters

n =9;
x= linspace(1,101,n+1)';
y = x;
z = x;

dx = x(2) - x(1);
dy = dx;
dz = dx;

x0 = x(1:n) + diff(x)/2;
y0 = y(1:n) + diff(y)/2;
%% Create data

% create cube of gravity equations
G = getG(x,y,z,x0,y0,dx,dy,dz);
% Play with anomaly here
m = zeros(n+1,n+1,n+1);
%m = m*exp(-((x-0.5*max(x) + (y - ))
m(1:2,8:9,3)=-8;
m(7:8,1:8,6)=7;
m(6,1:8,6)=5;
m(9,1:8,6)=5;
m(7:8,9,6) = 5;
m(3:8,1:4,8:9)=3;

% Random array generator
%mRandom = random('Normal',0,1,n+1,n+1,n+1);
load('mrandom.mat');
m = mRandom + m; %Turn addition of mRandom off to not include gaussian noise in data

% reshape M into vector
M = reshape(m,[length(G(:,1)),1]);
% Multiply to get D
D = G'*M;
% Reshape D into nxn grid for plotting
d = reshape(D,[n,n]);

%% Inversion

 % Call up W  - Currently GetW has 4 modes: "Eldad" mode, "Deep"  
 % mode, "Target" mode, and "Gradient" mode. 
 %Use getW(n+1,'mode','depth'); ('depth' not necessary if not
 %using target depth mode)
 depthTarget = 60;
 [W,Wz] = getW(n+1,depthTarget);
 
 % set W as a combination of given W and W for locating at depth
 WtW = W'*W;
 %WtW = W'*W;
 
 % Create an mref
mref = zeros(n+1,n+1,n+1);
% mref(2:9,3:7,4:5)=-1;
% mref(2:8,1:6,9)=2;
Mref = reshape(mref,[length(G(:,1)),1]);

% Choose an alpha if no alpha is selected
tol = 0.01; %      1 percent of sigmaN
[alphaDat, alphaValue] = alphaGen2(G,WtW,W,D,Mref,n,tol);
%alphaValue =1e6;

% Do inversion with nice alpha
MInv = (G*G' + alphaValue.*(WtW))\(G*D + alphaValue.*(WtW)*Mref);

% Reshape mInv
mInv = reshape(MInv,size(m));

%% PLOT

figure(1);
for ii = 1:n
 subP = subplot(3,3,ii);
        contourf(m(:,:,ii))
        title(sprintf('Depth z = %d',(ii-1)*10))
        caxis([min(min(min(m))) max(max(max(m)))])
        pos=get(subP, 'Position'); 
        set(subP, 'Position', [pos(1)+0.05 pos(2) pos(3) pos(4)]) 
        %cax = caxis; %Turn on for absolute color cross over between fig 1
        %and fig 2 (need to turn on corrisponding caxis in fig 2)
end
%rect = [left, bottom, width, height]

B=colorbar ('FontSize',12);
        set(B, 'Position', [.0314 .11 .0581 .8150]) 
         ylabel(B,'Gravity Anomaly [dimensionless]')
 
    
         
figure(2)
for ii = 1:n
subP = subplot(3,3,ii);
        contourf(mInv(:,:,ii))
        title(sprintf('Depth z = %d',(ii-1)*10))
        caxis([min(min(min(mInv))) max(max(max(mInv)))])
         pos=get(subP, 'Position'); 
        set(subP, 'Position', [pos(1)+0.05 pos(2) pos(3) pos(4)]) 
        %caxis(cax)  %Turn on for absolute color cross over between fig 1
        %and fig 2
end

B=colorbar ('FontSize',12);
         set(B, 'Position', [.0314 .11 .0581 .8150]) 
         ylabel(B,'Gravity Anomaly [dimensionless]')
   

figure(3)
contourf(d)
%B=colorbar ('FontSize',12);
       %  ylabel(B,'Gravity Anomaly [dimensionless]')
%set(gca,'YDir','reverse')