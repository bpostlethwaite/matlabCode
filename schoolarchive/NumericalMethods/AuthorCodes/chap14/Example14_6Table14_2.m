% Example 14.6 -- Table 14.2 -- differencing with uneven step

clear all
k = -5:1:-1;
h = .5*10.^k; h2 = 2*h;
exact = 1;

meth1 = h./(h.*h2)*exp(0) + (h./h2.*exp(h2) - h2./h.*exp(-h)) ./ (h+h2);
meth2 = (exp(h2)-exp(-h)) ./ (h+h2);

format long g
err1 = abs(meth1-exact)
err2 = abs(meth2-exact)