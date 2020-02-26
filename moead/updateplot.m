function updateplot(gen, pop, objs)

    str = sprintf('gen = %d', gen);

    hold off; 
    subplot(1, 2, 1);
    if size(objs, 2) == 2
        plot(objs(:, 1), objs(:, 2), 'ro', 'MarkerSize', 4);
        xlabel('f1', 'FontSize', 6);
        ylabel('f2', 'FontSize', 6);
    else
        plot3(objs(:, 1), objs(:, 2), objs(:, 3), 'ro', 'MarkerSize', 4);
        xlabel('f1', 'FontSize', 6);
        ylabel('f2', 'FontSize', 6);
        zlabel('f3', 'FontSize', 6);
    end
    title(str, 'FontSize', 8);
    box on;
    drawnow;

    subplot(1, 2, 2);
    if size(pop, 2) >= 3
        %plot3(pop(:, 1), pop(:, 2), pop(:, 3), 'ro', 'MarkerSize', 4);
        plot(pop(:, 1), pop(:, 2), 'ro', 'MarkerSize', 4);
        xlabel('x1', 'FontSize', 6);
        ylabel('x2', 'FontSize', 6);
        %zlabel('x3', 'FontSize', 6);    
    elseif size(pop, 2) >= 2
        plot(pop(:, 1), pop(:, 2), 'ro', 'MarkerSize', 4);
        xlabel('x1', 'FontSize', 6);
        ylabel('x2', 'FontSize', 6);
    end
    box on;
    drawnow;
    clear pareto objs pop;
end