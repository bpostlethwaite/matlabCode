%% basic Daisy World 
clear all
close all
clc

%% Set Variables

T_opt_b = 295; %Optimal growth temperature of black daisies in K (default = 295)
T_opt_w = 295; %Optimal growth temperature of white daisies in K (default = 295)
k_b = 17.5^(-2); %Growth rate bracketed between 35 K of T_opt
k_w = 17.5^(-2); %Growth rate bracketed between 35 K of T_opt
albedo_g = 0.5; %Albedo of bare ground
albedo_w = 0.75; %Albedo of white daisies
albedo_b = 0.25; %Albedo of black daisies
S = 1368/4; %Solar radiation in W/m2 received by a planet located 1 au from the Sun
q = 2.06*10^9; %Heat transfer coefficient to ensure thermal equilibrium (set so that q < 0.2*SL/sigma)
sigma = 5.67*10^-8;% Stefan Boltzman constant J/s/m2/K4
gamma = .3; %Death rate (default = 0.3)

%% Initial Values

ntime = 500; %Number of timesteps in main loop (default = 100)
Lset = [2 2.2 4.2]; %Luminosity (L is set to 2.5 where the black and white are 
%in equilibrium and similar values, 2.0 shows black rise, white fall)
hyst_max = Lset(3);
year = linspace(0,ntime,ntime);
temp = T_opt_w-49:1:T_opt_w+50;
alpha_w_store = zeros(1,ntime); %Storage vector for alpha_w
alpha_b_store = zeros(1,ntime); %Storage vector for alpha_b
temp_store = zeros(1,ntime); %Storage vector for Temperature
L_store = zeros(1,ntime); %Storage vector for luminosity
T_bare = zeros(1,ntime); % Creates a vector showing T without daisy influece
%% MAIN LOOP

%MAIN LOOPs%
% This first loop will cycle through 3 L conditions. The primary then 2
% others which are started at different L's. If the model did not show
% hysteresis then the values would match, but because it does display
% hysteresis the model depends on time, meaning the system has memory.

L = Lset;
for jj = 1:ntime
    T_bare(jj) = (S*L(1)*(1-albedo_g)/sigma)^(1/4);
    L_store(jj) = L(1);
    L(1) = L(1) + 0.01;   
end
    L_store=round(100*L_store);
    L_store = L_store/100;
for ii = 1:3
   L = Lset;
   T = T_bare(L_store==L(ii)); %Initial temperature of daisy world in K (default = 300)
   alpha_w = 0.01; %Fraction of land covered by white daisies
   alpha_b = 0.01; %Fraction of land covered by black daisies
   
for itime = 1:ntime
	alpha_g = 1 - alpha_b - alpha_w; %Fraction of land that is bare ground
    % The equation below set the main albedo as a function of hte
    % proportions of daisies coupled with their respective albedos
	A = alpha_w*albedo_w + alpha_b*albedo_b + albedo_g*(alpha_g); %mean planetary albedo
	%store variables for plotting%
    alpha_w_store(itime) = alpha_w; 
	alpha_b_store(itime) = alpha_b;
	temp_store(ii,itime) = T;
    % the equation below is a product of Stefan-Boltzmanns equation
    % S*L*(1-A) = sigma*T^4. In other words the temperature times sigma is
    % = to the solar energy actually being absorbed. It depends on the
    % current Albedo, which in turn is set by the proportions of Daisies.
	T = (S*L(ii)*(1-A)/sigma)^(1/4); %compute mean planetary temperature in radiative equilibrium
    
    % These next sets of equations create the temperature contribution from
    % each type of daisy and the bare planet. ie. Local temperatures.
	T_w = (q*(A - albedo_w) + T^4)^(1/4); %compute temperature of patch of white daisies
	T_b = (q*(A - albedo_b) + T^4)^(1/4); %compute temperature of patch of black daisies
	T_g = (q*(A - albedo_g) + T^4)^(1/4); %compute temperature of bare ground
    % The equations below set the proportions of daisies. These are based
    % on population replicator differential equations we have all seen in
    % our textbooks. ie. dw/dt=w*(p-w)*Birthrate-w*deathrate. In the code
    % form each change is iterated, so no differential equation is being
    % solved at once. so white daisies = white daisies + white daisies *
    % (proportion of land left*birthfunction(using local white daisy
    % temperature) - death rate)
	alpha_w = alpha_w + alpha_w*(alpha_g*betafn(T_w,T_opt_w,k_w) - gamma);
	alpha_b = alpha_b + alpha_b*(alpha_g*betafn(T_b,T_opt_b,k_b) - gamma);
    
    %% *#@! with the numbers here!
    %L=L*(1+sin(itime/4)/40);
    %L= L*(1+sin(pi*itime/40)/50);
    L(ii)=L(ii) + 0.01;
    
end    
end
%% Extra analysis
% Ok, so we need to shift the graphs of L(2) and L(3) over to so their
% starting point lines up with the appropriate Luminosity. The begining
% conditions, where they would otherwise be padded with zeros will be
% matched with T_bare.
% how much do we shift the graphs? where L_store == L(ii)
L = Lset;
temp_L1 = temp_store(1,:);
temp_L2 = [T_bare(1:find(L_store==L(2))-1),...
    temp_store(2,1:ntime+1-find(L_store==L(2)))];
temp_L3 = [T_bare(1:find(L_store==L(3))-1),...
    temp_store(3,1:ntime+1-find(L_store==L(3)))];


%% PLOTS
%DATA PLOTTING%

h=figure(1);
subplot(2,1,1)
    plot(L_store,temp_L1 - 273 ,'b',L_store,temp_L2 - 273,'r',...
        L_store,temp_L3 - 273,'k',L_store,T_bare -273,'--g');
    xlabel('Luminosity');
    ylabel('Temperature [C]');
    title('Release of daisies at staggered Luminosities')
    legend('Planet Temp w/ starting Luminosity','Temp w/ Luminosity II',...
        'Temp w/ Luminosity III', 'Barren Planet Temp','Location', 'Best');
    legend('boxoff');
    line([hyst_max,hyst_max ],...
        [temp_L1(L_store==hyst_max)-273, T_bare(L_store==hyst_max)-273],...
        'Color','k','Linestyle','--');
    text(hyst_max,temp_L1(L_store==hyst_max)-273,'\leftarrow Hysteresis loop',...
        'HorizontalAlignment','left','FontSize',12)
subplot(2,1,2)
    plot(year,L_store,'r');
    xlabel('year');
    ylabel('luminosity');
    title('Luminosity over time');
    

