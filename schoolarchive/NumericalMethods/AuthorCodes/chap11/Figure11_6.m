% Figure 11.6 : plot hat functions

clear all
close all
x = [0,1,2,4,5,7];

for i=1:7
  xx(i,:) = i-1:.01:i;
end
N = size(xx,2);

phi0 = zeros(size(xx));
phi0(1,:) = 1-xx(1,:);

xp = reshape(xx',1,N*7);
phi0 = reshape(phi0',1,N*7);

phi3 = zeros(size(xx));
phi3(3,:) = xx(3,:) - 2;
phi3(4,:) = 4 - xx(4,:);

phi3 = reshape(phi3',1,N*7);

phi4 = zeros(size(xx));
phi4(4,:) = xx(4,:) - 3;
phi4(5,:) = 5 - xx(5,:);

phi4 = reshape(phi4',1,N*7);
phi5 = zeros(size(xp));

plot(xp,phi0,xp,phi4,xp,phi3,xp,phi5)
hold on
plot(xp,phi3,'LineWidth',2,'Color',[.6,0,0])
     gtext('\phi_i')
     gtext('\phi_{i+1}')
     gtext('\phi_0')
     xlabel('x')
     ylabel('\phi')
