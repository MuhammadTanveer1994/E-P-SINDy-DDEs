function [bestTau, evalCount, elapsedtime, Theta_all, DX, X_data,bestpar3] = PSO(model, tau, par, DX_true, folderName)
 
tic;
Xi = [];
evalCount = 0;
paramHistory = [];
errorHistory = [];
par3History = [];  
iterCount = 0;

    function objective = optimizeParametersPSO(params)
        switch model
            case {'Rossler1', 'Rossler2'}
                tau(2) = params(1);   % Update tau_max
            case 'tau_3'
                tau(3) = params(1);  % Update tau_max
            case 'MG'
                tau = params(1);
                par(3) = params(2);
            otherwise
                tau = params(1);
        end

        [g, phi, par, tspan, X_data, DX, Xi, Theta_all, x_sindy, x_sol] = Popt(model, tau, par);

        err = norm(DX_true - Theta_all * Xi, 2) + norm(x_sol - x_sindy, 2);
        objective = err;

        evalCount = evalCount + 1;
        iterCount = iterCount + 1;

        paramHistory(iterCount, :) = params;
        errorHistory(iterCount) = err;

        switch model
            case 'MG'
                par3History(iterCount) = par(3);
        end
    end

switch model
    case {'Rossler1','Rossler2'}
        numDelays = 1;
        nVars = numDelays;
        lb = 1;
        ub = 2.5;
    case 'tau_3'
        numDelays = 1;
        nVars = numDelays;
        lb = [1];
        ub = [2.5];
    case 'MG'
        numDelays = 1;
        nVars = numDelays+1;
        lb = [0.1, 0.1];
        ub = [2, 20];
    otherwise
        numDelays = 1;
        nVars = numDelays;
        lb = [0.1];
        ub = [2];
end

options = optimoptions('particleswarm', ...
    'SwarmSize', 5, ...
    'HybridFcn', @patternsearch, ...
    'MaxIterations', 500, ...
    'Display', 'iter', ...
    'FunctionTolerance', 1e-3, ...
    'ObjectiveLimit', 1e-3, ...
    'PlotFcn', {@pswplotbestf});

[optParams, fval] = particleswarm(@optimizeParametersPSO, nVars, lb, ub, options);

switch model
    case 'MG'
        bestTau = optParams(1);
        bestpar3 = optParams(2);
        par(3) = bestpar3;
    otherwise
        bestTau = optParams(1);
end

fprintf('Total evaluateModel calls: %d\n', evalCount);
fprintf('Optimized tau: %f\n', bestTau);

figure;
switch model
    case 'MG'
        plot(paramHistory(:, 1), 'LineWidth', 2);
        hold on;
        plot(paramHistory(:, 2), 'LineWidth', 2);
        legend({'tau', 'par3'});
    otherwise
        plot(paramHistory(:, 1), 'LineWidth', 2);
        legend({'tau'});
end
title('Parameter Evolution');
xlabel('Iteration');
ylabel('Parameter Value');
grid on;
hold off;
savefig(fullfile(folderName, 'Parameter_Evolution.fig'));

figure;
plot(errorHistory, 'LineWidth', 2);
title('Error vs. Iteration');
xlabel('Iteration');
ylabel('Error');
grid on;
savefig(fullfile(folderName, 'Error_vs_Iteration.fig'));

figure;
switch model
    case 'MG'
        tauHistory = paramHistory(:,1);
        [minError, minIndex] = min(errorHistory);
        minTau = tauHistory(minIndex);
        minPar3 = par3History(minIndex);

        scatter(tauHistory, par3History, 50, errorHistory, 'filled');
        hold on;
        plot(minTau, minPar3, 'r*', 'MarkerSize', 15);
        xlabel('Tau');
        ylabel('par(3)');
        title('Parameter Space with Error (color)');
        colorbar;
        grid on;
        hold off;
        savefig(fullfile(folderName, 'MG_Parameter_Space.fig'));
    otherwise
        scatter(paramHistory(:, 1), errorHistory, 'filled');
        xlabel('Tau');
        ylabel('Error');
        title('Error vs Tau');
        grid on;
        savefig(fullfile(folderName, 'Error_vs_Tau.fig'));
end

elapsedtime = toc;

dataFileName = fullfile(folderName, 'PSO_Optimization_Data.mat');
save(dataFileName);
end
