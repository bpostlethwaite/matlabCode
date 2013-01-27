%% This script calls the function emissions.m and plots the CO2 emissions in giga-tons C per year.
%%This is the modified IPCC scenario A2 (see the pdf on the website)

%%clear memory and close figures
clear all
close all


%%time span in years
yr_MIN = 1800;
yr_MAX = 2150;

%% the function linspace makes a linearly spaced vector between the min and max years of n points.  n=200 below.
yr = linspace(yr_MIN, yr_MAX, 200);


%%Call the function emissions
e  = emissions(yr);



%%Plot the results
plot(yr, e);
xlabel('Time (yr)');
ylabel('CO_2 emissions (GtC/yr)');
title('Modified IPCC scenario A2')