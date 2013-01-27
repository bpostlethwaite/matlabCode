function beta = betafn_crittersC(ndaisies,t_a, T_opt, k)

if abs(t_a(1) - T_opt) < k^(-1/2)
    beta = 1 - k*((t_a - T_opt).^2);
    beta = [sum(beta(1:10)),sum(beta(11:20)),sum(beta(21:30))];
        
else
    beta = zeros(1,3);
end
