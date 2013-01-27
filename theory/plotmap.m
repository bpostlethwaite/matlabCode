function plotmap (obj, Ir, Ig, Ib)
    name = obj.type;
    rgb = obj.map;
    
    
    ax = linspace(0, 1, size(rgb, 1));
    
    plot(ax, rgb(:, 1), 'r', 'LineWidth',2.4)
    hold on
    plot(ax, rgb(:, 2), 'g', 'LineWidth',2.4)
    plot(ax, rgb(:, 3), 'b', 'LineWidth',2.4)
    plot(ax(Ir), rgb(Ir,1), 'ro','MarkerSize', 12, 'MarkerFaceColor', 'r')
    plot(ax(Ig), rgb(Ig,2), 'go','MarkerSize', 12, 'MarkerFaceColor', 'g')
    plot(ax(Ib), rgb(Ib,3), 'bo','MarkerSize', 12, 'MarkerFaceColor', 'b')
    title(sprintf('%s', name), 'FontSize', 16)
    grid on
    axis tight
    hold off
end