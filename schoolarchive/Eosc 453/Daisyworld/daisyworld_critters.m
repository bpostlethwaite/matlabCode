%% basic Daisy World 
clear all
close all
%% Set Variables 
% Regime Change == == set L to 2, and +0.001 to get dual regimes of reign
% control
% 
Lset = 2; %Luminosity (L is set to 2 where you get this cool spiral thing)
Y = Lset; %for renormalization later on
ntime = 2000; %Number of timesteps in main loop (default = 100)
ndaisies=30; %number of daisies in population, in numbers divisible by 3
T = 200; %Initial temperature of daisy world in K (default = 300)
T_opt_a = zeros(1,ndaisies); %create t_opt_a storage vector
T_opt_a(:) = 295; %Optimal growth temperature of daisies in K (default = 295)
k_a = 17.5^(-2); %Growth rate bracketed between 35 K of T_opt
k_B = 17.5^(-2); %Growth rate bracketing for Herbivores
albedo_a = 0:1/ndaisies:1-1/ndaisies; %Albedo of bare ground
albedo_g = 0.5; %Albedo of black daisies
S = 1368/4; %Solar radiation in W/m2 received by a planet located 1 au from the Sun
q = 2.06*10^9; %Heat transfer coefficient to ensure thermal equilibrium (set so that q < 0.2*SL/sigma)
sigma = 5.67*10^-8;% Stefan Boltzman constant J/s/m2/K4
gamma = .3; %Death rate (default = 0.3)
alpha_a = zeros(1,ndaisies); 
B = zeros(1,ndaisies/3);
C = zeros(1,3);
B_eaten = zeros(1,ndaisies);
C_eaten = zeros(1,ndaisies/3);
diversity = zeros(1,ntime);

%% Initial Values

alpha_a(:) = 0.01; %Fraction of land covered by daisies
B(:) = 0.01;
C(:) = 0.001;
year = linspace(0,ntime,ntime); %for plotting later on
temp = T_opt_a(1)-49:1:T_opt_a(1)+50; %for plotting later on
alpha_a_store = zeros(1,ndaisies); %Storage vector for daisies
B_store = zeros(1,ndaisies/3);
C_store = zeros(1,3);
temp_store = zeros(1,ntime); %Storage vector for Temperature
L_store = zeros(1,ntime); %Storage vector for luminosity
T_bare = zeros(1,ntime); % Creates a vector showing T without daisy influence

%% MAIN LOOP

L=Lset;

%MAIN LOOP%
for itime = 1:ntime
	alpha_g = 1 - sum(alpha_a); %Fraction of land that is bare ground
    alpha_Bg = 1 - sum(B);
    alpha_Cg = 1 - sum(C);
    % The equation below set the main albedo as a function of hte
    % proportions of daisies coupled with their respective albedos
	A = sum(alpha_a.*albedo_a) + albedo_g*(alpha_g); %mean planetary albedo
	%store variables for plotting%
    alpha_a_store = [alpha_a_store ; alpha_a]; 
    B_store = [ B_store ; B ];
    C_store = [ C_store ; C ];
	temp_store(itime) = T;
    % the equation below is a product of Stefan-Boltzmanns equation
    % S*L*(1-A) = sigma*T^4. In other words the temperature times sigma is
    % = to the solar energy actually being absorbed. It depends on the
    % current Albedo, which in turn is set by the proportions of Daisies.
	T = (S*L*(1-A)/sigma)^(1/4); %compute mean planetary temperature in radiative equilibrium
    T_bare(itime) = (S*L*(1-albedo_g)/sigma)^(1/4);
    % These next sets of equations create the temperature contribution from
    % each type of daisy and the bare planet. ie. Local temperatures.
	T_a = (q*(A - albedo_a) + T^4).^(1/4); %compute local temperature of daisies
	T_g = (q*(A - albedo_g) + T^4)^(1/4); %compute temperature of bare ground

 
    % The equations below set the proportions of daisies. These are based
    % on population replicator differential equations we have all seen in
    % our textbooks. ie. dw/dt=w*(p-w)*Birthrate-w*deathrate. In the code
    % form each change is iterated, so no differential equation is being
    % solved at once. 
	alpha_a = alpha_a + alpha_a.*(alpha_g*betafn_critters(ndaisies,T_a,T_opt_a,k_a) - gamma - B_eaten);
    %% HERBIVORES
    %turn daisies into food cells
    B_food = Bfood(alpha_a);
    B = B + B.*(alpha_Bg*B_food.*betafn_crittersB(ndaisies,T_a,T_opt_a,k_B) - gamma - C_eaten);
    
    % B_eaten is the B populations that eat certain daisies, basically
    % turning the smaller B vector into a representative larger alpha_a
    % vector
    B_eaten = beaten(B);
    
    %% CARNIVORES
    if itime == ntime*0.5
        C=C+0.1;
    else if itime > ntime*0.5
    C_food = [sum(B(1:3)),sum(B(4:7)),sum(B(8:10))];
    C = C + C.*(alpha_Cg*C_food.*betafn_crittersC(ndaisies,T_a,T_opt_a,k_B) - gamma);
    
    % B_eaten is the B populations that eat certain daisies, basically
    % turning the smaller B vector into a representative larger alpha_a
    % vector
    C_eaten = ceaten(C);
    
        end
    end
    
    
    %% *#@! with the numbers here!
% Note if you set the solar flare to +10 years and predator intro at one
% half of length you get two predators instead of usual one.
%     L_store(itime)= L; % luminosity over time
    
    if itime >= ntime*0.75 && itime <= ntime*0.75+10
        L = 4;
    else L=Lset;
    end

    
    %L= L*(1+sin(pi*itime/300)/1500);
    %L=L + 0.001;
%     if itime >= ntime-ntime/4-2 && itime <= ntime-ntime/4+2
%         L=4;
%     else L = L*(1+sin(pi*itime/40)/50);
%     end
    
    
    
   
end
%% Extra analysis

%     %Diversity Index - - - Calculated via Shannon Index.
%     % ie. Sum of the all the proportions of individuals in species a,b,c... etc to total species
%     % times the natural log of that proportion.
%     alpharound=(round(10*alpha_a_store))/10;
%     %diversity(itime)=-sum((alpharound/sum(alpharound)).*log(alpharound/sum(alpharound)));
%     for ii = 1:ntime
%     
%         H = (alpharound(ii,1)/sum(alpharound(ii,1)))*log(alpharound(ii,1)/sum(alpharound(ii,1)));
%     for jj=2:ndaisies
%         Y =(alpharound(ii,jj)/sum(alpharound(ii,jj)))*log(alpharound(ii,jj)/sum(alpharound(ii,jj)));
%         H=H+Y;
%     end
%         diversity(ii)=-H;
%     end

%% PLOTS
%DATA PLOTTING%
figure(1);
subplot(3,1,1)
plot(L_store,temp_store,L_store,T_bare,'--');
%plot(year,temp_store,year,T_bare);
xlabel('Luminosity');
ylabel('Temperature [K]');
legend('Temp w/ daisies','Temp w/o daisies','Location', 'Best');
legend('boxoff');
%xlim([1,ntime]);

subplot(3,1,2)
xlim([1 ntime])
hold on
plot(alpha_a_store,'g')
plot(B_store,'b')
plot(C_store,'r')
hold off
ylabel('fractional cover')
xlim([1,ntime]);
legend('daisies','herbivores','carnivores')

subplot(3,1,3);
plot(year,L_store,'r',year,diversity,'b');
xlabel('year');
ylabel('luminosity');

figure(2);
subplot(3,1,1)
plot(alpha_a_store);

subplot(3,1,2)

plot(B_store);
legend('herb1','herb2','herb3','herb4','herb5','herb6','herb7','herb8','herb9','herb10')
subplot(3,1,3)
plot(C_store);
legend('Carnivore 1','Carnivore 2', 'Carnivore 3');
