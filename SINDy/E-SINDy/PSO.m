function [bestTau, evalCount, elapsedtime, Theta_all, DX,bestpar3] = PSO(model, tau, par, DX_true, folderName)

    tic;
    Xi = [];
    evalCount = 0;
    paramHistory = [];  
    errorHistory = []; 
    tau1History = [];  
    tau2History = [];   
    tau3History = [];  
    tauHistory = []; 
    par3History = [];  
    iterCount = 0;     
    
    switch model
        case {'Rossler1','Rossler2'}
            bestTau1 = [];
            bestTau2 = [];
            bestTau = [bestTau1, bestTau2]; 
           
        case 'tau_3'
            bestTau1 = [];
            bestTau2 = [];
            bestTau3 = [];
            bestTau = [bestTau1, bestTau2, bestTau3]; 
        case 'MG'
            bestTau = [];
            bestpar3 = [];
        otherwise
            bestTau = [];
    end

    function objective = optimizeParametersPSO(params)
        switch model
            case {'Rossler1', 'Rossler2'}
                tau1 = params(1);
                tau2 = params(2);
                tau = [tau1, tau2];
     
            case 'tau_3'
                tau1 = params(1);
                tau2 = params(2);
                tau3 = params(3);
                tau = [tau1, tau2, tau3];
            case 'MG'
                tau = params(1);
                par(3) = params(2); 
            otherwise
                tau = params(1);
        end

        [g, phi, par, tspan, X_data, DX, Xi, Theta_all, x_sindy, x_sol] = Eopt(model, tau, par);

        err = norm(DX_true - Theta_all * Xi, 2) + norm(x_sol - x_sindy, 2);
        objective = err;

        evalCount = evalCount + 1;

        iterCount = iterCount + 1;
        paramHistory(iterCount, :) = params;
        errorHistory(iterCount) = err;

        switch model
            case {'Rossler1', 'Rossler2'}
                tau1History(iterCount) = tau1;
                tau2History(iterCount) = tau2;
            case 'tau_3'
                tau1History(iterCount) = tau1;
                tau2History(iterCount) = tau2;
                tau3History(iterCount) = tau3;
            case 'MG'
                tauHistory(iterCount) = tau;
                par3History(iterCount) = par(3);
            otherwise
                tauHistory(iterCount) = tau;
        end
    end

    switch model
        case 'Rossler2'
            numDelays = 2;
            nVars = numDelays;  
            lb = [0.1, 1];    
            ub = [1.5, 3];         
            
        case 'Rossler1'
           
            numDelays = 2;
            nVars = numDelays;  
            lb = [0.000001, 1];    
            ub = [1, 3];        
        case 'tau_3'
            numDelays = 3;
            nVars = numDelays;  
            lb = [1, 1, 1];    
            ub = [2, 2.5, 2.5];         
        
        case 'MG'
            numDelays = 1;
            nVars = numDelays + 1;  % One for tau and one for par(3)
            lb = [0.1, 0.1];      
            ub = [2, 20];         
        otherwise
            numDelays = 1;
            nVars = numDelays;      
            lb = 0.1;            
            ub = 1.5;              
    end

    options = optimoptions('particleswarm', ...
        'SwarmSize', 5, ...
        'HybridFcn', @patternsearch, ...
        'MaxIterations', 1000, ...
        'Display', 'iter', ...
        'FunctionTolerance', 1e-3, ...
        'ObjectiveLimit', 1e-3, ...
        'PlotFcn', @pswplotbestf);

    [optParams, fval] = particleswarm(@optimizeParametersPSO, nVars, lb, ub, options);

    switch model
        case {'Rossler1', 'Rossler2', 'TDN'}
            bestTau1 = optParams(1);
            bestTau2 = optParams(2);
            bestTau = [bestTau1, bestTau2];
        case 'tau_3'
            bestTau1 = optParams(1);
            bestTau2 = optParams(2);
            bestTau3 = optParams(3);
            bestTau = [bestTau1, bestTau2, bestTau3];
        case 'MG'
            bestTau = optParams(1);
            bestpar3 = optParams(2);
            par(3) = bestpar3;  
        otherwise
            bestTau = optParams(1);
    end

    fprintf('Total evaluateModel calls: %d\n', evalCount);

    switch model
        case {'Rossler1', 'Rossler2'}
            fprintf('Best Tau1: %f, Best Tau2: %f\n', bestTau1, bestTau2);
        case 'tau_3'
            fprintf('Best Tau1: %f, Best Tau2: %f, Best Tau3: %f\n', bestTau1, bestTau2, bestTau3);   
        case 'MG'
            fprintf('Best Tau: %f, Best par3: %f\n', bestTau, bestpar3);
        otherwise
            fprintf('Best Tau: %f\n', bestTau);
    end

    % Common figures for all models
    figure;
    plot(errorHistory, 'LineWidth', 2);
    title('Error vs. Iteration');
    xlabel('Iteration');
    ylabel('Error');
    grid on;
    savefig(fullfile(folderName, 'Error_vs_Iteration.fig'));

    figure;
    hold on;
    paramLabels = {};
    switch model
        case {'Rossler1', 'Rossler2'}
            plot(paramHistory(:, 1), 'LineWidth', 2);  
            plot(paramHistory(:, 2), 'LineWidth', 2);  
            paramLabels = {'tau1', 'tau2'};
        case 'MG'
            plot(paramHistory(:, 1), 'LineWidth', 2);  
            plot(paramHistory(:, 2), 'LineWidth', 2);  
            paramLabels = {'tau', 'par3'};
        case 'tau_3'
            plot(paramHistory(:, 1), 'LineWidth', 2); 
            plot(paramHistory(:, 2), 'LineWidth', 2);  
            plot(paramHistory(:, 3), 'LineWidth', 2);  
            paramLabels = {'tau1', 'tau2', 'tau3'};
        otherwise
            plot(paramHistory(:, 1), 'LineWidth', 2); 
            paramLabels = {'tau'};
    end
    title('Parameter Evolution');
    xlabel('Iteration');
    ylabel('Parameter Value');
    legend(paramLabels);
    grid on;
    hold off;
    savefig(fullfile(folderName, 'Parameter_Evolution.fig'));

    % Model-specific figures
    if numDelays == 2
        [minError, minIndex] = min(errorHistory);
        minTau1 = tau1History(minIndex);
        minTau2 = tau2History(minIndex);
        
         figure;
    scatter(tau1History, errorHistory, 'filled');
    xlabel('Tau1');
    ylabel('Error');
    title('Error vs Tau1');
    grid on;
    savefig(fullfile(folderName, 'Error_vs_Tau1.fig'));

    figure;
    scatter(tau2History, errorHistory, 'filled');
    xlabel('Tau2');
    ylabel('Error');
    title('Error vs Tau2');
    grid on;
    savefig(fullfile(folderName, 'Error_vs_Tau2.fig'));

        figure;
        scatter3(tau1History, tau2History, errorHistory, 15, errorHistory, 'filled');
        hold on
        plot3(minTau1, minTau2, minError, 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r'); 
        xlabel('Tau1');
        ylabel('Tau2');
        zlabel('Error');
        title('Error Surface vs Tau1 and Tau2');
        colorbar;
        grid on;
        savefig(fullfile(folderName, 'Error_Surface_Tau1_Tau2.fig'));

         figure;
    scatter3(tau1History, tau2History, errorHistory, 15, errorHistory, 'filled');
    hold on;
    plot3(minTau1, minTau2, minError, 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
    xlabel('Tau1');
    ylabel('Tau2');
    title('Error Surface vs Tau1 and Tau2 with Minimum Error Point');
    colorbar;
    view(2);   
    hold off;
    savefig(fullfile(folderName, 'Error_Surface_Tau1_Tau2_TopDown.fig'));


    elseif numDelays == 3
        [minError, minIndex] = min(errorHistory);
        minTau1 = tau1History(minIndex);
        minTau2 = tau2History(minIndex);
        minTau3 = tau3History(minIndex);   

        % Main comprehensive figure for 3-delay case
        figure;
        subplot(2,2,1);
        scatter(tau1History, tau2History, 50, errorHistory, 'filled');
        xlabel('\tau_1');
        ylabel('\tau_2');
        title('\tau_1 vs \tau_2 (Error as color)');
        colorbar;
        hold on;
        plot(tau1History(minIndex), tau2History(minIndex), 'r*', 'MarkerSize', 15);
        hold off;
        
        subplot(2,2,2);
        scatter(tau2History, tau3History, 50, errorHistory, 'filled');
        xlabel('\tau_2');
        ylabel('\tau_3');
        title('\tau_2 vs \tau_3 (Error as color)');
        colorbar;
        hold on;
        plot(tau2History(minIndex), tau3History(minIndex), 'r*', 'MarkerSize', 15);
        hold off;
        
        subplot(2,2,3);
        scatter(tau1History, tau3History, 50, errorHistory, 'filled');
        xlabel('\tau_1');
        ylabel('\tau_3');
        title('\tau_1 vs \tau_3 (Error as color)');
        colorbar;
        hold on;
        plot(tau1History(minIndex), tau3History(minIndex), 'r*', 'MarkerSize', 15);
        hold off;
        
        subplot(2,2,4);
        scatter3(tau1History, tau2History, tau3History, 50, errorHistory, 'filled');
        hold on;
        plot3(minTau1, minTau2, minTau3, 'r*', 'MarkerSize', 15);
        xlabel('\tau_1');
        ylabel('\tau_2');
        zlabel('\tau_3');
        title('3D Parameter Space (Error as color)');
        colorbar;
        view(45, 30);
        hold off;
        
        sgtitle('Multi-view Analysis of Tau Parameters and Error');
        savefig(fullfile(folderName, 'Tau_MultiView_Analysis.fig'));
        
        figure;
        scatter3(tau1History, tau2History, tau3History, 50, errorHistory, 'filled');
        hold on;
        plot3(tau1History(minIndex), tau2History(minIndex), tau3History(minIndex), ...
            'r*', 'MarkerSize', 15, 'DisplayName', 'Minimum Error Point');
        xlabel('\tau_1');
        ylabel('\tau_2');
        zlabel('\tau_3');
        title('Tau Recovery Space with Error (color)');
        colorbar;
        legend('Parameter Space', 'Minimum Error Point');
        grid on;
        view(45, 30);
        
        text(tau1History(minIndex), tau2History(minIndex), tau3History(minIndex), ...
            sprintf('\nMinimum Error Point\n\\tau_1 = %.2f\n\\tau_2 = %.2f\n\\tau_3 = %.2f\nError = %.2e', ...
            tau1History(minIndex), tau2History(minIndex), tau3History(minIndex), minError), ...
            'FontSize', 10);
            
        savefig(fullfile(folderName, 'Tau3_Error_Surface.fig'));


    elseif strcmp(model, 'MG')
        [minError, minIndex] = min(errorHistory);
        minTau = tauHistory(minIndex);
        minPar3 = par3History(minIndex);
        
        figure;
        scatter(tauHistory, par3History, 50, errorHistory, 'filled');
        hold on;
        plot(minTau, minPar3, 'r*', 'MarkerSize', 15);
        xlabel('Tau');
        ylabel('par(3)');
        title('Parameter Space with Error (color)');
        colorbar;
        grid on;
        savefig(fullfile(folderName, 'MG_Parameter_Space.fig'));
    else
        figure;
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
