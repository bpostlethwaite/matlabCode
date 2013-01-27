% Example 15.13 : using the adaptive routine with different tol values

[Qe,mesh,fevals] = quads(0,1,1.e-13);     % ``exact value''

exps = 3:2:11;
tols = 10.^(-exps);
format long g

for j = 1:length(tols)
  [Q,mesh,fevals] = quads(0,1,tols(j));
  tabvals = [abs(Qe-Q),fevals]
end
    

