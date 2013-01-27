clear all
close all
load PCA.mat

y = [x1;x2;x3;x4];
ymat = y';
t = 1:40;

% Mean reduce data
yb = y - repmat(mean(y,2),1,size(y,2));

% Covariance Matrix
C = 1/(39)*(yb*yb');

% Eigs and Eigfunctions of Covariance 
[V,D]= eig(C);
[D,order] = sort(diag(D),'descend');  %# sort eigenvalues in descending order
V = V(:,order);

% Principal components, temporal coefficients amplitudes... (diff names 
% for a,
% Obtained by projecting mean reduced data vector yb onto eigenvector ej
a = V'*yb;
dotsum = sum(a.*a,2);
varpercent = dotsum./sum(dotsum)*100;
totalVarY = sum(dotsum)*1/length(x1);

% Grab Power spectrum of 'a'
for ii = 1:4;
[Pxx(ii,:),F(ii,:)] = periodogram(a(ii,:),[],[],1);
end


% Deprojected principal components back into time series, using first two
% modes (first two a and eigs)
Y2 = V(:,1:2)*a(1:2,:);
Y3 = V(:,1:3)*a(1:3,:);

figure(1)
for ii = 1:4
    subplot(4,2,2*ii-1)
        plot(t,a(ii,:))
        legend('principal component','Location','NorthWest')
        legend boxoff
        title(sprintf('Time Series for PC Mode %d \n contributing %2.1f%% variance in y',ii,varpercent(ii)))
        xlabel('t')
        ylabel('x')
    subplot(4,2,2*ii)
        plot(F,Pxx(ii,:))
        title(sprintf('Mode %i Periodogram / Power spectrum',ii))
        xlabel('dimensionless frequency')
        ylabel('Amplitude')
end

figure(2)
for ii = 1:4
    subplot(4,1,ii)
        plot(t,yb(ii,:),t,Y2(ii,:),'r.',t,Y3(ii,:),'gs')
        legend('Time Series','Modes 1 & 2','Modes 1,2 & 3',...
            'Location','BestOutside')
        legend boxoff
        title(sprintf(['Time Series for X%i \n Points show degree'...
            'of variation captured by increasing modes'],ii))
end


