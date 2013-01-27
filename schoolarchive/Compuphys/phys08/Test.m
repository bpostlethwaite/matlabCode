clear all
close all

% N = 300;
% 
% 
% %a = 1.58;
% aa = linspace(1,1.9,50);
% for jj = 1:length(aa)
%    a = aa(jj); 
%    sum1 = 0;
% sig1 = 0;
% sum2 = 0;
% sig2 = 0;
% 
% %a = 1.586
% for ii = 1:N
%     lm    = -a*exp(-1) + a;
%     lm2   = -a*exp(-pi) + a;
%     x     = rand(1) * lm;
%     x2    = rand(1) * lm2;
%     yx    = -log(abs(1-x/a));
%     yx2   = -log(abs(1-x2/a));
%     fx    = exp(-yx^2)/(a*exp(-yx)); 
%     fx2   = 1/(a*exp(-yx2)*(yx2^2 + cos(yx2)^2));
%     sum1  = sum1 + fx;
%     sum2  = sum2 + fx2;
%     sig1  = sig1 + fx^2;
%     sig2  = sig2 + fx2^2;
% end
% 
% Monte(jj)   = sum1*lm/N;
% Monte2(jj)  = sum2*lm2/N;
% SIGma       = sig1*lm/N;
% SIGma2      = sig2*lm2/N;
% VAR(jj)     = abs(SIGma - Monte(jj)^2);
% VAR2(jj)    = abs(SIGma2 - Monte2(jj)^2);
% end



A      = load('phys08.dat');
aa     = A(:,1);
Monte  = A(:,2);
Monte2 = A(:,3);
VAR    = A(:,4);
VAR2   = A(:,5);


subplot(2,1,1)
    plot(aa,VAR,'*',aa,Monte,'+')
    legend('abs|Variance|','Integral','Location','Best')
    title('For exp(-x^2) - Variance and Integral')
    xlabel('value of ''a''')
    ylabel('function value')
subplot(2,1,2)    
    plot(aa,VAR2,'*',aa,Monte2,'+')
    legend('abs|Variance|','Integral','Location','Best')
    title('For 1/(x^2 + cos(x)^2) - Variance and Integral')
    xlabel('value of ''a''')
    ylabel('function value')

