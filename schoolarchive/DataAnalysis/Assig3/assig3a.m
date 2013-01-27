.
% Data Analysis
clear all
close all

load spec.mat
X = [x1;x2];
%% Question 1

figure(1)
for ii = 1:2
[P(ii,:),F(ii,:)] = periodogram(X(ii,:),[],[],1);

m(1,:) =  [11.42, 6.88, 1.122];
m(2,:) =  [11.63, 2.85, 2.45];

varpercent = m(ii,:)/sum(P(ii,:))*100;

maxfreq = F(ii,P(ii,:)==m(ii,1));
    subplot(2,2,2*ii-1)
        plot(t,X(ii,:))
        title(sprintf('time series for x%1.0i',ii))
        xlim([0,t(end)])
        ylabel('amplitude')
        xlabel('time in Days')
    subplot(2,2,2*ii)
        plot(F(ii,:),P(ii,:))
        title(sprintf(['AutoSpectrum for x%1.0i with 3 peak frequency showing'...
            '%2.1f%% %2.1f%% %2.1f%% variation'],ii,varpercent(1),varpercent(2),varpercent(3)))
        xlim([0,0.15])
        ylabel('amplitude')
        xlabel('frequency in Days')
        
end
    
%% Question 2
lag = 51;
n = length(x1);

for ii = 1:2
    for jj = 1:lag
        L(jj,:) = X(ii,jj:n-lag+jj);
    end
    C = 1/length(L(:,1)) * L*L';
    [V,D] = eig(C);
    [D,order] = sort(diag(D),'descend');
    fprintf('First 6 RC''s in time series %i capture:  \n',ii)
    fprintf('RC %i: %2.1f %% of the variation\n',[[1:6];[100*D(1:6)./sum(D)]'])
    V = V(:,order);
    
    a = V'*L;
    
    
    
    DELAY = 1;
    for qq = 1:6
        Y(:,:,qq) = V(:,qq)*a(qq,:);
        [Pxx(qq,:),F(qq,:)] = periodogram(Y(DELAY,:,qq),[],[],1);

    end
     
    Tlength = length(Y(1,:,1));
    tt = t(1:Tlength);

    figure(5*ii)
    ind = 1:10:length(Y(:,1,1));
    for k = 1:length(ind)
        legendmatrix{k}=strcat('lag series: ',num2str(ind(k)));
    end
    plot(tt,Y(ind,:,1))
    legend(legendmatrix)
    title(sprintf(['Approximated data from principal RC of Time Series %i\n' ...
        'Plotting multiple Lagged Series'],ii))

    
    figure(3*ii)
    for jj = 1:3;
        RC = jj;
        subplot(3,2,2*jj-1)
            plot(tt,Y(DELAY,:,2*RC-1),':',tt,Y(DELAY,:,2*RC),'--')
            title(sprintf('Time Series %i approximation with mode %i and %i',ii,2*RC-1,2*RC))
            xlim([0,tt(end)])
            ylabel('Amplitude')
            xlabel('Time in Days')
            legend(sprintf('mode %i',2*RC-1),sprintf('mode %i',2*RC),'Location','Best')
        subplot(3,2,2*jj)
            plot(F(jj,:),Pxx(2*RC-1,:),':',F(jj,:),Pxx(2*RC,:),'--')
            title(sprintf('Power Spectrum of modes %i and %i',ii,2*RC-1,2*RC))
            legend(sprintf('mode %i',2*RC-1),sprintf('mode %i',2*RC),'Location','Best')
            xlim([0,0.15])
            ylabel('amplitude')
            xlabel('frequency in Days')
    end

     
    figure(4)
        subplot(2,1,ii)
        plot(V(:,1))
        title(sprintf('Principal SSA Eigenvector for time series %i',ii))
        xlim([0,length(V(:,1))])
        xlabel('Component of EigenVector')
        ylabel('Magnitude of Component')

end
