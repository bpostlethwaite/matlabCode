% Example 2.2 -- Figure 2.2

t = 0:.002:1;
tt = exp(-t) .* (sin(2*pi*t)+2);
rt = single(tt);
round_err = (tt - rt) ./tt ;
plot (t,round_err,'b-');
title('error in sampling exp(-t)(sin(2\pi t)+2) in single precision')
xlabel('t')
ylabel('roundoff error')

% relative error should be about eta = eps(single)/2 
rel_round_err = max(abs(round_err)) / (eps('single')/2)