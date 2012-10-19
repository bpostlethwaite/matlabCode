% Example 14.1 -- Table 14.1 -- Figure 14.1 : simple difference methods

clear all
clf
k = -8:1:-1;
h = 10.^k; h2 = 2*h;
exact = 1;

meth1 = (exp(h)-1) ./ h;
meth2 = (exp(h)-exp(-h)) ./ h2;
meth3 = (-exp(h2)+8*exp(h)-8*exp(-h)+exp(-h2)) ./ (12*h);

format long g
err1 = abs(meth1-exact)
err2 = abs(meth2-exact)
err3 = abs(meth3-exact)

loglog(h,err1,h,err2,h,err3)
xlabel('h')
ylabel('absolute error')
legend('1st order','2nd order','4th order')