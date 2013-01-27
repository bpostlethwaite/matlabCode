function [roots,iter] = rootgrabber(func,nprobe,a,b,tol,maxiter,print)
% rootgrabber finds roots of the given function.
% First probes function nprobe times for zero crossings on interval [a,b]
% and then performs root finding via the secant method and returns result
% within tolerance tol, doing a max of maxiter iterations. If print = true,
% the last option, it also prints out roots + iteration.

% Right now the secant method employs a computationally expensive trick to
% get at the boundaries. It calculates duplicates right now, but throws
% them away. This needs to change in upcoming version.





try
[roots,iter] = secant(func,a,b,nprobe,tol,maxiter);
[roots, m, n] = unique(roots);
iter = iter(m);
catch ME    
end

if print == true;
    for ii= 1:length(roots)
        fprintf('found root: %3.6g with %2d iterations\n',roots(ii),iter(ii))
    end
end


