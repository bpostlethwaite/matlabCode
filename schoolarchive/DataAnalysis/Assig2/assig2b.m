clear all
close all
load CCA.mat

X = [x1;x2;x3]';
Y = [y1;y2;y3]';

[A,B,r,U,V] = canoncorr(X,Y);

F = cov(X)*A;
G = cov(Y)*B;

% Projections of 1st, then 1st two, then all 3 modes.
for ii = 1:3
XX(:,:,ii) = U(:,1:ii)*F(:,1:ii)';
YY(:,:,ii) = V(:,1:ii)*G(:,1:ii)';
end

t = 1:size(X,1);

% SVD (max covariance)
Cxy = 1/(size(X,1)-1)*(X'*Y);
[U2,S2,V2] = svd(Cxy);

gph = char('b','r','g');


% GET Power Spectrum
for ii = 1:3
[Puu(ii,:),Fu(ii,:)] = periodogram(U(:,ii),[],[],1);
[Pvv(ii,:),Fv(ii,:)] = periodogram(V(:,ii),[],[],1);
end



figure(1)
for ii = 1:3
    hold on
    figure(1)
    subplot(3,1,ii)
        plot(t,X(:,ii),t,Y(:,ii),'r')
        legend(sprintf('x%i',ii),sprintf('y%i',ii))
        title(sprintf('Time Series for x%i and y%i',ii,ii))
    figure(2)
        subplot(3,2,2*ii-1)
            plot(t,U(:, ii),t,V(:,ii),'r')
            title(sprintf('canonical variates u%i and v%i',ii,ii))
            legend(sprintf('u%i',ii),sprintf('v%1',ii))
        subplot(3,2,2*ii)
            plot(Fv(ii,:),Puu(ii,:),Fv,Pvv(ii,:),'r')
            title(sprintf('Power Spectrum for u%i & v%i',ii,ii))
            xlabel('dimensionless freq')
            ylabel('Amplitude')
            legend(sprintf('u%i',ii),sprintf('v%i',ii))
end
  

    figure(12)
for ii=1:3

    hold on
    subplot(3,1,ii)
        plot(t,X(:,ii),t,XX(:,ii,1),'r*',t,XX(:,ii,2),'g*',t,XX(:,ii,3),'k*')
legend('original Tseries','inverse mapping using 1st vector',...
    'inverse mapping using 1st & 2nd vector','inverse mapping all vectors',...
    'Location','NorthEastOutside')
end


Xv=[A(1,1:2),U2(1,1:2)];
Yv=[A(2,1:2),U2(2,1:2)];
Zv=[A(3,1:2),U2(3,1:2)];
z = zeros(1,4);
Gph = char(gph,'k');

figure(123)
for ii = 1:4
    hold on
    quiver3(0,0,0,Xv(ii),Yv(ii),Zv(ii),Gph(ii))
    grid on
end
    
    
    