function beta = betafn_crittersB(ndaisies,t_a, T_opt, k)

if abs(t_a(1) - T_opt) < k^(-1/2)
    beta = 1 - k*((t_a - T_opt).^2);
    beta = [sum(beta(1:3)),sum(beta(4:6)),sum(beta(7:9)),sum(beta(10:12)),sum(beta(13:15)),...
        sum(beta(16:18)),sum(beta(19:21)),sum(beta(22:24)),sum(beta(25:27)),sum(beta(28:30))];
else
    beta = zeros(1,ndaisies/3);
end
