function e = emissions(yr)

%%Function defining A2 emission scenario over the interval 1990-2100
%%extended to pre-industrial (assuming linear increase from 0 in 1850 to
%%1900) and assuming full cessation of CO_2 input at 2101

%%For additional information see http://www.grida.no/climate/ipcc/emission

t_yr =      [0 1850 1990  2000   2010  2020 2030   2040  2050 2060  2070  2080   2090  2100 2110 2120 10000];

%%CO2 forcing stops at 2100: note that need some zeros close in time to
%%approximate stepwise shutoff with linear interp

%%
e_GtC_yr = [0  0   6.875 8.125 9.375  12.5 14.375 16.25 17.5 19.75 21.25 23.125 26.25 28.75 0     0     0];




e = interp1(t_yr , e_GtC_yr , yr , 'linear');