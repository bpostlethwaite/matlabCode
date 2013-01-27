function beta = betafn(temp, T_opt, k)
if abs(temp - T_opt) < k^(-1/2)
    beta = 1 - k*(temp - T_opt)^2;
else
    beta = 0;
end
