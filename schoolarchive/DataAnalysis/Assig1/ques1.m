clear all
close all

load corr.mat;

X = [x,x2,x3];
Y = [y,y2,y3];

TYPE = {'Pearson';'Spearman'};
for i = 1:3
    Xi = [ones(length(X(:,i)),1),X(:,i)]; % Xi = 1 + xi
    Yi = Y(:,i); % Yi = yi
    Beta = (Xi'*Xi)\Xi'*Yi; %Solve for coefficents Beta (a0 a1)
    R(:,i) = Xi*Beta;       % Regression Line = Xi*Beta
    for j = 1:2
        rho(i,j) = corr(X(:,i),Y(:,i),'type',TYPE{j}); % Get Correlation
        varNR(i,j) = (1-rho(i,j)^2);
    end
    
    % New Regression line using slope of spearman corrleation coefficients
    
    R2(:,i) = Xi*[Beta(1);rho(i,2)*(std(Y(:,i))/std(X(:,i)))];
    
    sprintf('Regression Slope is: %1.2f\nCorrelation Coefficient * std dev ratio is: %1.2f\n',...
        Beta(2),rho(i,1).*(std(Y(:,i))/std(X(:,i))))
    
end    


% Spearman data sorter
Ys = sort(Y);
Xs = sort(X);

    for jj = 1:length(Ys(1,:))
        for ii = 1:length(Ys(:,1))
            ind = find(Ys(ii,jj)==Y(:,jj));
            Ysp(ii,jj) = ind/(length(Ys(:,1))-1);
        end
    end
    
    for jj = 1:length(Xs(1,:))
        for ii = 1:length(Xs(:,1))
            ind = find(Xs(ii,jj)==X(:,jj));
            Xsp(ii,jj) = ind;
        end
    end
    
for ii = 1:3    
figure(23)
    subplot(3,1,ii)
    plot(X(:,ii),Ysp(:,ii),'*r',X(:,ii),Y(:,ii)./max(Y(:,ii)),'sb')
end


TITLE = {'x & y scatterplot', 'x2 & y2 scatterplot','x3 & y3 scatterplot'};
for i = 1:3
    figure(1)
    h=subplot(3,1,i);
        plot(X(:,i),Y(:,i),'r*',X(:,i),R(:,i),'b')
        title(TITLE(i))
        hleg = legend('x at y','regression','Location','EastOutside');
        legend boxoff
        text(0.5,0.25,sprintf('%s Correlation: %1.2f', TYPE{1},rho(i,1)),...
            'Units','normalized')
        text(0.5,0.1,sprintf('%s Correlation: %1.2f',TYPE{2}, rho(i,2)),...
            'Units','normalized')
        set(h,'TickDir','out')
figure(2)        
    h=subplot(3,1,i);
        plot(X(:,i),Y(:,i),'r*',X(:,i),R2(:,i),'b')
        title(TITLE(i))
        hleg = legend('x at y','regression','Location','EastOutside');
        legend boxoff
        text(0.5,0.25,sprintf('%s Correlation: %1.2f', TYPE{1},rho(i,1)),...
            'Units','normalized')
        text(0.5,0.1,sprintf('%s Correlation: %1.2f',TYPE{2}, rho(i,2)),...
            'Units','normalized')
        set(h,'TickDir','out')
end
%%
%Reguratate above calcs but set the spearmanafied X and Y values
%to X and Y and repeat calcs

% X = Xsp;
% Y = Ysp;
% 
% for i = 1:3
%     Xi = [ones(length(X(:,i)),1),X(:,i)]; % Xi = 1 + xi
%     Yi = Y(:,i); % Yi = yi
%     Beta = (Xi'*Xi)\Xi'*Yi; %Solve for coefficents Beta (a0 a1)
%     R(:,i) = Xi*Beta;       % Regression Line = Xi*Beta
%     for j = 1:2
%         rho(i,j) = corr(X(:,i),Y(:,i),'type',TYPE{j}); % Get Correlation
%         varNR(i,j) = (1-rho(i,j)^2);
%     end
%     
%     % New Regression line using slope of spearman corrleation coefficients
%     
%     R2(:,i) = Xi*[Beta(1);rho(i,2)];
%     
%     sprintf('Regression Slope is: %1.2f\nCorrelation Coefficient * std dev ratio is: %1.2f\n',...
%         Beta(2),rho(i,1).*(std(Y(:,i))/std(X(:,i))))
%     
% end    
% 
% TITLE = {'x & y scatterplot', 'x2 & y2 scatterplot','x3 & y3 scatterplot'};
% for i = 1:3
%     figure(5)
%     h=subplot(3,1,i);
%         plot(X(:,i),Y(:,i),'r*',X(:,i),R(:,i),'b')
%         title(TITLE(i))
%         hleg = legend('x at y','regression','Location','EastOutside');
%         legend boxoff
%         text(0.5,0.25,sprintf('%s Correlation: %1.2f', TYPE{1},rho(i,1)),...
%             'Units','normalized')
%         text(0.5,0.1,sprintf('%s Correlation: %1.2f',TYPE{2}, rho(i,2)),...
%             'Units','normalized')
%         set(h,'TickDir','out')
% figure(6)        
%     h=subplot(3,1,i);
%         plot(X(:,i),Y(:,i),'r*',X(:,i),R2(:,i),'b')
%         title(TITLE(i))
%         hleg = legend('x at y','regression','Location','EastOutside');
%         legend boxoff
%         text(0.5,0.25,sprintf('%s Correlation: %1.2f', TYPE{1},rho(i,1)),...
%             'Units','normalized')
%         text(0.5,0.1,sprintf('%s Correlation: %1.2f',TYPE{2}, rho(i,2)),...
%             'Units','normalized')
%         set(h,'TickDir','out')
% end
