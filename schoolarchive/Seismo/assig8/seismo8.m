%Lab8
close all
clear all
run quakes;

for p = 1:2
    vel=6;
    x=linspace(0,100,100);
    y=linspace(0,100,100);
    [XX,YY] = meshgrid(x,y);
    
    E = zeros(100,100);
    T = E;
    K = E;
    tq{1} = tq1;
    tq{2} = tq2;
    
    % Compute epsilon
    for ii = 1:length(tq{p})
        T = T + (tq{p}(ii)-(sqrt((XX-xcrd(ii)).^2+(YY-ycrd(ii)).^2))/vel);
    end
    T = T/length(tq{p});
    for ii = 1:length(tq{p})
        E = E + (tq{p}(ii)-(sqrt((XX-xcrd(ii)).^2+(YY-ycrd(ii)).^2)/vel)-T).^2;
    end
    % Find min points
    [Y,X] = find(E==min(min(E)));
    
    %Compute sigma square
    Tbest  = mean(tq{p}-(((X-xcrd).^2+(Y-ycrd).^2).^(1/2))/vel);
    sigma  = 1/10*sum((tq{p} -(sqrt((X-xcrd).^2+(Y-ycrd).^2)/vel)-Tbest).^2);
    
    %Chi squared
    for ii=1:length(xcrd)
        K = K + (tq{p}(ii)-(sqrt((XX-xcrd(ii)).^2+(YY-ycrd(ii)).^2)/vel)-T).^2;
    end
    K = K ./ sigma;
    % Plot
    figure (p)
    imagesc(E)
    hold on
    contour (K,[18.31 18.31],'k')
    plot(xcrd,ycrd,'o','MarkerSize',10,...
        'MarkerFaceColor','g','MarkerEdgeColor','k')
    plot(X,Y,'o','MarkerSize',10,...
        'MarkerFaceColor','r','MarkerEdgeColor','k')
    hold off
end