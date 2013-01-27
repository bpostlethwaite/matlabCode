clear all 
close all


load NLR.mat

meanRanF = mean(ytestranf(:,41:80),2);
meanBays = mean(ytestbays,2);

ind = 1:80;

MIN = min([min(Bayrms),min(RanForestrms)]);
MAX = max([max(Bayrms),max(RanForestrms)]);

for ii = 1:80
    C(ii) = corr(ytestbays(:,ii),ytestranf(:,ii));
    V(ii) = var(ytestbays(:,ii)-ytestranf(:,ii));
end

%
for ii = 1:80
    figure(345)
    subplot(2,1,1)
        plot3(x1test,x2test,ytestbays(:,ii),'ro',...
              x1test,x2test,ytestranf(:,ii),'b^')
          xlim([-1,1])
          ylim([-1,1])
          zlim([-2.5 2.5])
          grid on
          axis square
          
    subplot(2,1,2)
        plot(ind,Bayrms,ind,RanForestrms,ind,C)
        line([ii ii],[MIN MAX],'LineWidth',4,'Color',[.4 .9 .8])
        legend('Baysian','Random Forest','Correlation')
    pause(1)
end
%}
%%

%{
 plot3(x1test,x2test,meanBays,'ro',...
              x1test,x2test,meanRanF,'b^')
          xlim([-1,1])
          ylim([-1,1])
          zlim([-2.5 2.5])
          grid on
          axis square
%}

%csvwrite('BenPostle_ytests.csv',M)
% Rany 42
% Bays 16