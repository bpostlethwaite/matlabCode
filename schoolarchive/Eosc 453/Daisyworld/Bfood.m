function out = Bfood(alpha_a)
% BFOOD - turns white daisies into tasty veggie treats (it contracts
% daisies into a 1/3 population as to be = to herbavore vector size
out = [sum(alpha_a(1:3)),sum(alpha_a(4:6)),sum(alpha_a(7:9)),sum(alpha_a(10:12)),sum(alpha_a(13:15)),...
        sum(alpha_a(16:18)),sum(alpha_a(19:21)),sum(alpha_a(22:24)),sum(alpha_a(25:27)),sum(alpha_a(28:30))];
end
