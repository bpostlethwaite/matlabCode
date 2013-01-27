% test
clear all, close all, clc

v = [4 6 8 10 12 15 20]; % velocities
u = 1./v; % slownesses
p = linspace(.01,.999,500)*u(1); % p values for all take off angles
lp = length(p);
dz = 3;
% r(:,1) = ones(1,lp);

% calculate Xp
X = zeros(lp,1);
T = X;
R = X;

for i = 1:length(v)-1
    r(:,i) = (p < u(i));
    R = R + r(:,i);
    eta = (u(i)^2 - p.^2).^.5;
    Xp(:,i)  = 2*dz*r(:,i)'.*(p./eta);
    Tp(:,i)  = r(:,i).*((2*dz*u(i)^2)./eta');
    X = X + Xp(:,i);
    T = T + Tp(:,i);
end
%% 
xm = [0 max(X)];
ym = [0 max(T)];
figure(1), hold on
plot(linspace(0,100),linspace(0,100)*u(1))
for i = 1:lp-1
    plot(X(end-i:end),T(end-i:end)) 
    xlim(xm), ylim(ym)
    
    pause(.0005)
    
end
hold off
%%
t = round(T*5);
x2 = X/2;

for i = 1:lp
rayx{i} = linspace(0,X(i)/2,t(i));
rayz{i} = -linspace(0,R(i)*dz,t(i));
end

for i = 1:lp
rayz2{i} = -rayz{i} - R(i)*dz;
rayx2{i} = x2(i) + rayx{i};
end

figure(2), hold on
xlim([0 40]), ylim([-12.1 0.1])
for i = 1:lp
plot(rayx{i},rayz{i})
plot(rayx2{i},rayz2{i})
end




