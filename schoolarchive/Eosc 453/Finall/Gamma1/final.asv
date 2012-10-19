

E = 1e9;
te = linspace(100,20000,1000);
v = 0.25;
rr = linspace(0,10000,1000);
r = 3000;
phi = 0.4;
% Uses the right te2 for the 0 to 5km lid thickness
te2=[1,1000,2000,3000,4000,5000,10000]*phi;
g = 1.314;
rho = 917;
delta_rho = 4;
lamda = 4*r;
lamda2 = 4*rr;
lamda3 = 4*linspace(1000,5000,5);
x=linspace(0,20000,1000);
xx = linspace(-2000,2000,1000);
c3=10;
Vo=2*r*delta_rho*g;



D = E*te.^3/(12*(1-v^2));
h = 2*r*delta_rho/rho*1./(1+16*D*pi^4./(lamda^4*rho*g));
lid=te2/phi/1000;


for ii = 1:length(te2)

D2(ii) = E*te2(ii)^3/(12*(1-v^2));
h2(ii,:) = 2*rr.*delta_rho/rho*1./(1+16*D2(ii)*pi^4./(lamda2.^4*rho*g));
alpha(ii) = (4*D2(ii)/(rho*g)).^(1/4);
w(ii,:)= 10*exp(-x./alpha(ii)).*(cos(x./alpha(ii))+sin(x./alpha(ii)));

end

%% calculate viscosity

gs= 2e-5;
p = 2;
q = 59e3;
RR = 8.314;
A = 42*(1.97e-5)*(9.1e-4)/(250*8.314);
c = (gs^p)/(A);
T = 128;

visc = c/(1) * exp(q/(RR*T));

relax = 2*visc/E

diapirforce = 2*delta_rho*g/(r*pi);

visc2 = c/(1) * exp(q/(RR*250));
risetime = visc2/diapirforce
%%
%Calculate total ice depth with ratios

Depth= te2/Rt;

figure(1)
plot(te/1000,h)

figure(2)
plot(rr/1000,h2)
ylim([0,50])
xlim([0,10])
xlabel('Diapir radius [km]')
ylabel('surface deformation amplitude [m]')
legend('0 km','1 km','2 km','3 km','4km','5km','10 km','Location','Best')
title('surface deformation versus diapir radius for increasing stagnant lid thickness')
line([0 9],[10,10],'LineWidth',2,'color',[.8 .8 .8])


figure(3)
plot(x/1000,w)
title('Amplitude versus lenticulae radius')
legend('0 km','1 km','2 km','3 km','4km','5km','10 km','Location','Best')
ylabel('surface deformation amplitude [m]')
xlabel('lenticulae radius')
ylim([-2,10])
xlim([0,20])


% if the minimum detectable dome amplitude is about 10m, then the
% minimum radius is about 1 to 2 times this. 