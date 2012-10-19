% Figure 11.7 : plot Hermite basis functions
%
clear all
close all

x = [0,1,2,4,5,7];

for i=1:7;
  xx(i,:) = i-1:.01:i;
end
N = size(xx,2);

xi0 = zeros(size(xx));
xi0(1,:) = 1 - 3*xx(1,:).^2 + 2*xx(1,:).^3;

xp = reshape(xx',1,N*7);
xi0 = reshape(xi0',1,N*7);

xi3 = zeros(size(xx));
xi3(4,:) = 1 - 3*(xx(4,:)-3).^2 + 2*(xx(4,:)-3).^3;
xi3(3,:) = 3*(xx(3,:)-2).^2 - 2*(xx(3,:)-2).^3;

xi3 = reshape(xi3',1,N*7);

xi4 = zeros(size(xx));
xi4(5,:) = 1 - 3*(xx(5,:)-4).^2 + 2*(xx(5,:)-4).^3;
xi4(4,:) = 3*(xx(4,:)-3).^2 - 2*(xx(4,:)-3).^3;

xi4 = reshape(xi4',1,N*7);
xi5 = zeros(size(xp));

subplot(2,1,1)
plot(xp,xi0,xp,xi4,xp,xi3,xp,xi5)
hold on
plot(xp,xi3,'LineWidth',2,'Color',[.6,0,0])
     gtext('\xi_k')
     gtext('\xi_{k+1}')
     gtext('\xi_0')
     xlabel('x')
     ylabel('\xi')
     
%%%%%%%%%%%%%%%%%%%

ei0 = zeros(size(xx));
ei0(1,:) = xx(1,:) - 2*xx(1,:).^2 + xx(1,:).^3;

xp = reshape(xx',1,N*7);
ei0 = reshape(ei0',1,N*7);

ei3 = zeros(size(xx));
ei3(4,:) = xx(4,:)-3 - 2*(xx(4,:)-3).^2 + (xx(4,:)-3).^3;
ei3(3,:) = -(xx(3,:)-2).^2 + (xx(3,:)-2).^3;

ei3 = reshape(ei3',1,N*7);

ei4 = zeros(size(xx));
ei4(5,:) = xx(5,:)-4 - 2*(xx(5,:)-4).^2 + (xx(5,:)-4).^3;
ei4(4,:) = -(xx(4,:)-3).^2 +(xx(4,:)-3).^3;

ei4 = reshape(ei4',1,N*7);
ei5 = zeros(size(xp));

subplot(2,1,2)
plot(xp,ei0,xp,ei4,xp,ei3,xp,ei5)
hold on
plot(xp,ei3,'LineWidth',2,'Color',[.6,0,0])
     gtext('\eta_k')
     gtext('\eta_{k+1}')
     gtext('\eta_0')
     xlabel('x')
     ylabel('\eta')
     
