function hat = hatfunc(T,T_opt,k)

if 1 - (2*abs(T - T_opt))/k == abs(T-T_opt)*k/2
    hat = 1 - (2*abs(T - T_opt))/k;
else
    hat = 0;
end