function plotResults(sol, solSINDy, t_r, t_ae, x_tn, T, N, n_a, model, dist, folderName)

    figure;
     t_all = linspace(0, T, 1000);
    o_sol = deval(sol, t_all);
    i_sindy = deval(solSINDy, t_all);

    splitIndex = round(length(t_all) * t_r);
    splitTime = t_all(splitIndex); 

    plot(t_all, o_sol, 'Color', [0.7 0.7 0.7], 'LineWidth', 3);
    hold on;
    plot(t_all, i_sindy, 'k--', 'LineWidth', 2);
     plot(t_ae(1:length(x_tn)), x_tn, 'ro', 'LineWidth', 2); 
    yLimits = ylim;
    line([splitTime splitTime], yLimits, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 1.5); % Vertical line
    hold off;

    xlabel('t', 'FontSize', 13, 'Interpreter', 'latex');
    set(gca,'fontsize',20,'fontname','times');
    figName = fullfile(folderName, sprintf('SimulationResults-%s.fig', model));
    savefig(figName);
    if strcmp(model, 'Rossler2')
        figure; 
        plot3(sol.y(1,:), sol.y(2,:), sol.y(3,:), 'Color', [0.7 0.7 0.7], 'LineWidth', 3);
        hold on;
        plot3(solSINDy.y(1,:), solSINDy.y(2,:), solSINDy.y(3,:), 'k--', 'LineWidth', 2);
        hold off;
        xlabel('x', 'FontSize', 13, 'Interpreter', 'latex');
        ylabel('y', 'FontSize', 13, 'Interpreter', 'latex');
        zlabel('z', 'FontSize', 13, 'Interpreter', 'latex');

        set(gca, 'fontsize', 20, 'fontname', 'times');
    figName = fullfile(folderName, sprintf('SimulationResults-3D-%s.fig', model));
    savefig(figName);
    elseif strcmp(model, 'Rossler1')
        figure; 
        plot3(sol.y(1,:), sol.y(2,:), sol.y(3,:), 'Color', [0.7 0.7 0.7], 'LineWidth', 3);
        hold on;
        plot3(solSINDy.y(1,:), solSINDy.y(2,:), solSINDy.y(3,:), 'k--', 'LineWidth', 2);
        hold off;
        xlabel('x', 'FontSize', 13, 'Interpreter', 'latex');
        ylabel('y', 'FontSize', 13, 'Interpreter', 'latex');
        zlabel('z', 'FontSize', 13, 'Interpreter', 'latex');

        set(gca, 'fontsize', 20, 'fontname', 'times');
    figName = fullfile(folderName, sprintf('SimulationResults-3D-%s.fig', model));
    savefig(figName);
    elseif strcmp(model, 'tau_3')
    figure; 
    plot(sol.y(1,:), sol.y(2,:), 'Color', [0.7 0.7 0.7], 'LineWidth', 3);
    hold on;
    plot(solSINDy.y(1,:), solSINDy.y(2,:), 'k--', 'LineWidth', 2);
    hold off;
    xlabel('y_1', 'FontSize', 13, 'Interpreter', 'latex');
    ylabel('y_2', 'FontSize', 13, 'Interpreter', 'latex');
    set(gca, 'fontsize', 20, 'fontname', 'times');
    figName = fullfile(folderName, sprintf('PhasePlot-%s.fig', model));
    savefig(figName);
    end


end

