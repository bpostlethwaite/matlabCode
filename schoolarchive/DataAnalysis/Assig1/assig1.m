clear all
close all

load corr.mat;

X = [x,x2,x3];
Y = [y,y2,y3];

TYPE = {'Pearson';'Spearman'};
for i = 1:3
    Xi = [ones(length(X(:,i)),1),X(:,i)]; % Xi = 1 + xi
    Yi = Y(:,i); % Yi = yi
<<<<<<< HEAD
    Beta = (Xi'*Xi)\(Xi'*Yi); %Solve for coefficents Beta (a0 a1)
=======
    Beta = (Xi'*Xi)\Xi'*Yi; %Solve for coefficents Beta (a0 a1)
>>>>>>> origin/master
    R(:,i) = Xi*Beta;       % Regression Line = Xi*Beta
    for j = 1:2
        rho(i,j) = corr(X(:,i),Y(:,i),'type',TYPE{j}); % Get Correlation
    end
    
    sprintf('Regression Slope is: %1.2f\nCorrelation Coefficient * std dev ratio is: %1.2f\n',...
        Beta(2),rho(i,1).*(std(Y(:,i))/std(X(:,i))))
end    

TITLE = {'x & y scatterplot', 'x2 & y2 scatterplot','x3 & y3 scatterplot'};
for i = 1:3
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
end



