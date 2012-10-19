%Lab8 Ben Postlethwaite
%close all
%clear all
run quakes;

for p = 1:2
    vel=6;
    n = 101;
    x=linspace(0,100,n);
    y=linspace(0,100,n);
    [XX,YY] = meshgrid(x,y);
    E = zeros(n,n);
    T = E;
    K = E;
    E2 = E;
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
    sigma  = 1/10*sum((tq{p} -(sqrt((XX(Y,X)-xcrd).^2+(YY(Y,X)-ycrd).^2)/vel)-T(Y,X)).^2);
    %Chi squared
    for ii=1:length(xcrd)
        K = K + (tq{p}(ii)-(sqrt((XX-xcrd(ii)).^2+(YY-ycrd(ii)).^2)/vel)-T).^2;
        check = (tq{p}(ii)-(sqrt((XX-xcrd(ii)).^2+(YY-ycrd(ii)).^2)/vel)-T).^2;
    end
    K = K ./ sigma;
    
    % Plot
    figure (p+a)
    imagesc(y,x,E)
    hold on
    contour (x,y,K,[18.31 18.31],'k','LineWidth',3)
    plot(xcrd,ycrd,'o','MarkerSize',10,...
        'MarkerFaceColor','g','MarkerEdgeColor','k')
    plot(XX(Y,X),YY(Y,X),'o','MarkerSize',10,...
        'MarkerFaceColor','r','MarkerEdgeColor','k')
    hold off
    fprintf('For Earthquake number %d\n',p)
    fprintf('Best Location is: %2.2d x and %2.2d y\n',XX(Y,X),YY(Y,X))
    fprintf('Best Origin Time is: %2.2f s\n',T(Y,X))
    fprintf('"Variance" epsilon at this location is: %2.2f\n',E(Y,X))
    fprintf('sigma^2 at this location is: %2.2f\n',sigma)
    fprintf('chi^2 at this location is: %2.2f\n',K(Y,X))
    fprintf('\n\n')
end

