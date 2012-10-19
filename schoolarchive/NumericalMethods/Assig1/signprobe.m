function subi = signprobe(func,a,b,nprobe)
% SUBI probes func on subinterval and returns a 2xn matrix of subintervals
% Probes nprobe times looking for change of sign. Returns NaN if no
% interval found.

subint = a:b/nprobe:b;
fx = func(subint);
subi=[NaN;NaN];
jj=1;

for ii = 1:length(subint)-1
    if fx(ii)*fx(ii+1) <= eps
        subi(:,jj) = [subint(ii);subint(ii+1)];
        jj = jj + 1;
    end
end
