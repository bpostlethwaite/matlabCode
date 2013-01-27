clear all
close all
h = 40000;
rho1 = 2700;
rho2 = 3300;
beta1 = 3500;
beta2 = 4500;
%c = [3800,4000,4200,4400];
c = linspace(beta1,beta2,1000);
p = 1./c;
mu1 = rho1*beta1^2;
mu2 = rho2*beta2^2;
w = linspace(0,2,1000);

for jj = 1:length(p)
    f1(jj,:) = tan(h*w.*sqrt(1/beta1^2 - p(jj)^2));
    f2(jj,1) = (mu2 * sqrt(p(jj)^2 - 1/beta2^2))/(mu1 * sqrt(1/beta1^2 - p(jj)^2));
    W{jj} = getequal(f1(jj,:),f2(jj),w);
end


if length(c) < 6
    fprintf('\n\nC value         Period\n')
    for ii = 1:4
        fprintf('%d m          %1.2f s\n',c(ii),2*pi/W{ii}(1))
        F2(ii,:) = ones(1,length(w)).*f2(ii);
    end

figure(1)
  plot(w,f1,'*')
    hold on
  plot(w,F2)
    ylim([-2,4])
    legend([repmat('C = ',4,1) num2str(c')])

    
else
   % plot(W{1:10}(1))
    
end


