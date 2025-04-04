function plotResult(solA, solSINDy, t_r,t_ae,tspan, x_tn,T,N,M,n_a,model,dist,n,folderName)
    figure;

    t_all = linspace(tspan(1),tspan(end),1000);
    o_sol = deval(solA, t_all);
    i_sindy = deval(solSINDy, solSINDy.x);

    splitIndex = round(length(t_all) * t_r);
    splitTime = t_all(splitIndex); 

    plot(t_all, o_sol(1:n,:), 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5); 
    hold on;
    plot(solSINDy.x, i_sindy(1:n,:), 'k--', 'LineWidth', 2); 

    yLimits = get(gca, 'ylim'); 
    line([splitTime splitTime], yLimits, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 1.5); % Vertical line

    hold on;
    t_ae = t_ae(1:length(x_tn)); 
    plot(t_ae, x_tn(1:n,:), 'ro', 'LineWidth', 2);

    hold off;

    xlabel('t', 'FontSize', 13, 'Interpreter', 'latex');
    distLeg = ['Dist=' num2str(dist)];
    set(gca,'fontsize',20,'fontname','times');
    hold off;

    figName = fullfile(folderName, sprintf('SimulationResults-%s.fig', num2str(model)));
savefig(figName); 

if ismember(model, {'Rossler1','Rossler2'})
        figure; 
        plot3(solA.y(1,:), solA.y(2,:), solA.y(3,:), 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
        hold on;
        plot3(solSINDy.y(1,:), solSINDy.y(2,:), solSINDy.y(3,:), 'k--', 'LineWidth', 2);
        hold off;
        xlabel('x', 'FontSize', 13, 'Interpreter', 'latex');
        ylabel('y', 'FontSize', 13, 'Interpreter', 'latex');
        zlabel('z', 'FontSize', 13, 'Interpreter', 'latex');

        set(gca, 'fontsize', 20, 'fontname', 'times');
    figName = fullfile(folderName, sprintf('SimulationResults-3D-%s.fig', model));
    savefig(figName);
    end
end
