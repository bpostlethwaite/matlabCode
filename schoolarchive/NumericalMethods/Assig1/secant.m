function [roots,iter] = secant(f,a,b,nprobe,tol,maxiter)
% SECANT uses the secant method to find the roots of function f
% Finds along interval a,b probed initially for changes in sign nprobe
% times and roots are found to within tolerance tol. Number of iterations
% is also returned

subi = signprobe(f,a,b,nprobe);

if isnan(subi(1)) || isnan(subi(2)) 
    fprintf('No Sign Change (root) along specified interval\n')
end

for jj = 1:length(subi)
    x = subi(:,jj);
    itr = 1;
    while (abs(f(x(end)))) > tol && (abs(x(end)-x(end-1)) > tol*(1+abs(x(end))));
        x(end+1) = x(end) - f(x(end))*(x(end)-x(end-1)) / (f(x(end))-f(x(end-1)));
        if (abs(f(x(end))) >= 0.5*abs(f(x(end-1)))) % || (iter > maxiter)
            x = signprobe(f,subi(1,jj),subi(2,jj),3); % If change in x is small,
            % ie. if we are in a 'flat' interval of function, reduce the
            % interval around the sign change and restart Newtons method
        end
        if (itr > maxiter) % Break out of while loop if we breach max iter
            x = NaN;
            break
        end
        itr = itr + 1;
    end
    roots(jj) = x(end);
    iter(jj) = itr;
end
