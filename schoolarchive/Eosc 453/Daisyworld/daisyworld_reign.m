% Daisyworld Rein Control


%%  Variables
albedo_w = 1;
albedo_b = 0;
albedo_bed = 0.5;
gamma = 0.0001;
T_opt = 295;
S = 1368/4;
L = 1;

%% Initial Values
ntime = 2000;
alpha_w = 0.01;
alpha_b = 0.01;
T_w = 280;
T_b = 280;
D = 0.5;
k = 17.5^(-2);
% alpha_w_store = zeros(1:ntime);
% alpha_b_store = zeros(1:ntime);
% T_w_store = zeros(1:ntime);
% T_b_store = zeros(1:ntime);


%% Loop
for itime = 1:ntime

alpha_w_store(itime) = alpha_w;
alpha_b_store(itime) = alpha_b;
T_w_store(itime) = T_w;
T_b_store(itime) = T_b;

A_w = albedo_w*alpha_w + (1-alpha_w)*albedo_bed;
A_b = albedo_b*alpha_w + (1-alpha_b)*albedo_bed;

T_w = S*L*(1-A_w)-D*(T_b-T_w);
T_b = S*L*(1-A_b)-D*(T_w-T_b);

d_alpha_w = hatfunc(T_w,T_opt,k);
d_alpha_b = hatfunc(T_b,T_opt,k);

alpha_w = (1-gamma)*alpha_w + gamma*d_alpha_w;
alpha_b = (1-gamma)*alpha_b + gamma*d_alpha_b;


end