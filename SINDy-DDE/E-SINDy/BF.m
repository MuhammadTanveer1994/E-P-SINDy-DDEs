function [bestTau, evalCount, elapsedtime, Theta_all, DX, X_data,bestPar3] = BF(model, tau, DX_true, par, folderName)
 
    if nargin < 4 || isempty(par)
        par = [];  
    end
    if nargin < 5 || isempty(folderName)
        folderName = pwd;  
    end

    tic;
    evalCount  = 0;
    Xi         = [];
    Theta_all  = [];
    DX         = [];
    X_data     = [];

    switch model
        case 'Rossler1'
            tau1Grid = linspace(0.000001,1, 100);  
            tau2Grid = linspace(0.1, 4, 100);  

        case 'Rossler2'
            tau1Grid = linspace(0.1, 1.5, 100);   
            tau2Grid = linspace(1, 3, 100); 

        case 'MG'
            tauGrid  = linspace(0.1, 2, 10);   
            par3Grid = linspace(0.1, 20, 10);  

        otherwise
            tauGrid = linspace(0.1, 1.5, 100);
    end


    bestError  = inf;
    bestParams = struct('tau', [], 'par3', []);
    
    switch model
        case {'Rossler1', 'Rossler2'}

            errors2D = zeros(length(tau1Grid), length(tau2Grid));

            for i = 1:length(tau1Grid)
                for j = 1:length(tau2Grid)
                    [g, phi, par, tspan, X_data, DX, Xi, Theta_all] = ...
                        Eopt(model, [tau1Grid(i), tau2Grid(j)]);
                    evalCount = evalCount + 1;
                    currentError = norm(DX_true - Theta_all * Xi, 2);%+ norm(x_sol   - x_sindy,       2)
                    errors2D(i,j) = currentError;
                    if currentError < bestError
                        bestError      = currentError;
                        bestParams.tau = [tau1Grid(i), tau2Grid(j)];
                    end
                end
            end

            bestTau = bestParams.tau;
            
            fprintf('\n*** Optimal Parameters for %s ***\n', model);
            fprintf('  Tau1 = %.4f, Tau2 = %.4f\n', bestTau(1), bestTau(2));
            fprintf('  Best Error = %.6e\n', bestError);
            fprintf('  Total Eopt calls = %d\n', evalCount);

            figure('Name', sprintf('%s Parameter Search', model));
            [Tau1Mesh, Tau2Mesh] = meshgrid(tau1Grid, tau2Grid);

            contourf(Tau1Mesh, Tau2Mesh, errors2D', 20, 'LineStyle', 'none');
            hold on;
            plot(bestTau(1), bestTau(2), 'r*', 'MarkerSize', 10, 'LineWidth', 2);
            colorbar;
            xlabel('\tau_1');
            ylabel('\tau_2');
            title(sprintf('Error Distribution for %s', model));
            grid on;

            figName = fullfile(folderName, sprintf('error-%s.fig', model));
            savefig(figName);

        case 'MG'
            errors2D = zeros(length(tauGrid), length(par3Grid));

            paramMG = zeros(length(tauGrid)*length(par3Grid), 3);
            idx = 1;

            for i = 1:length(tauGrid)
                for j = 1:length(par3Grid)
                    [g, phi, par, tspan, X_data, DX, Xi, Theta_all] = ...
                        Eopt(model, tauGrid(i), par3Grid(j));

                    evalCount = evalCount + 1;

                    currentError = norm(DX_true - Theta_all * Xi, 2);%+ norm(x_sol   - x_sindy,       2)

                    errors2D(i,j) = currentError;

                    paramMG(idx,1) = tauGrid(i);
                    paramMG(idx,2) = par3Grid(j);
                    paramMG(idx,3) = currentError;
                    idx = idx + 1;

                    if currentError < bestError
                        bestError       = currentError;
                        bestParams.tau  = tauGrid(i);
                        bestParams.par3 = par3Grid(j);
                    end
                end
            end

            bestTau = bestParams.tau;
            bestPar3 = bestParams.par3;
            fprintf('\n*** Optimal Parameters for MG ***\n');
            fprintf('  Tau  = %.4f\n', bestTau);
            fprintf('  par3 = %.4f\n', bestParams.par3);
            fprintf('  Best Error = %.6e\n', bestError);
            fprintf('  Total Eopt calls = %d\n', evalCount);

            figMG = figure;
            [TauMesh, Par3Mesh] = meshgrid(tauGrid, par3Grid);

            subplot(1,2,1);
            contourf(TauMesh, Par3Mesh, errors2D', 20, 'LineStyle', 'none');
            hold on;
            plot(bestParams.tau, bestParams.par3, 'r*', 'MarkerSize', 12, 'LineWidth', 2);
            colorbar;
            xlabel('\tau');
            ylabel('par3');
            title('Error Contour & Parameter Distribution');
            grid on;

            subplot(1,2,2);
            surf(TauMesh, Par3Mesh, errors2D', 'EdgeColor','none');
            hold on;
            plot3(bestParams.tau, bestParams.par3, bestError, 'r*', 'MarkerSize', 12, 'LineWidth', 2);
            colorbar;
            xlabel('\tau');
            ylabel('par3');
            zlabel('Error');
            title('Error Surface');
            grid on;
            view(45, 30);

            figName = fullfile(folderName, 'MG_optimization_landscape.fig');
            savefig(figName);

            [minErrorVal, minIndex] = min(paramMG(:,3));
            minTau  = paramMG(minIndex,1);
            minPar3 = paramMG(minIndex,2);

            figure('Name','MG: Tau vs par3 (Error color)');
            scatter(paramMG(:,1), paramMG(:,2), 50, paramMG(:,3), 'filled');
            hold on;
            plot(minTau, minPar3, 'r*', 'MarkerSize', 12, 'LineWidth', 2);
            xlabel('\tau');
            ylabel('par3');
            title('MG Parameter Space (\tau, par3) with Error as Color');
            c = colorbar;
            c.Label.String = 'Error';
            grid on;
            savefig(fullfile(folderName, 'MG_Tau_par3_Scatter.fig'));

            figure;
            scatter3(paramMG(:,1), paramMG(:,2), paramMG(:,3), 50, paramMG(:,3), 'filled');
            hold on;
            plot3(minTau, minPar3, minErrorVal, 'r*', 'MarkerSize', 12, 'LineWidth', 2);
            xlabel('\tau');
            ylabel('par3');
            zlabel('Error');
            title('MG 3D Parameter Space (\tau, par3, error)');
            c = colorbar;
            c.Label.String = 'Error';
            grid on;
            view(45, 30);
            savefig(fullfile(folderName, 'MG_Tau_par3_Error_3Dscatter.fig'));

            figure;
            subplot(1,2,1);
            scatter(paramMG(:,1), paramMG(:,3), 40, paramMG(:,3), 'filled');
            hold on;
            plot(minTau, minErrorVal, 'r*', 'MarkerSize', 12, 'LineWidth', 2);
            xlabel('\tau');
            ylabel('Error');
            title('Error vs \tau');
            colorbar; grid on;

            subplot(1,2,2);
            scatter(paramMG(:,2), paramMG(:,3), 40, paramMG(:,3), 'filled');
            hold on;
            plot(minPar3, minErrorVal, 'r*', 'MarkerSize', 12, 'LineWidth', 2);
            xlabel('par3');
            ylabel('Error');
            title('Error vs par3');
            colorbar; grid on;

            savefig(fullfile(folderName, 'MG_Error_vs_Tau_par3.fig'));

        
        otherwise
            errors = zeros(length(tauGrid), 1);

            for i = 1:length(tauGrid)
                [g, phi, par, tspan, X_data, DX, Xi, Theta_all] = ...
                    Eopt(model, tauGrid(i));

                evalCount = evalCount + 1;

                currentError = norm(DX_true - Theta_all * Xi, 2);%+ norm(x_sol   - x_sindy,       2)

                errors(i) = currentError;

                if currentError < bestError
                    bestError      = currentError;
                    bestParams.tau = tauGrid(i);
                end
            end

            bestTau = bestParams.tau;
            fprintf('\n*** Optimal Parameter for %s ***\n', model);
            fprintf('  Tau = %.4f\n', bestTau);
            fprintf('  Best Error = %.6e\n', bestError);
            fprintf('  Total Eopt calls = %d\n', evalCount);

            figure('Name', sprintf('%s Parameter Search', model));
            plot(tauGrid, errors, 'o-');
            hold on;
            plot(bestTau, bestError, 'r*', 'MarkerSize', 10, 'LineWidth', 2);
            xlabel('\tau');
            ylabel('Error');
            title(sprintf('Error vs. \\tau for %s', model));
            grid on;

            figName = fullfile(folderName, sprintf('error-%s.fig', model));
            savefig(figName);
    end

   bestTau = bestParams.tau;
bestPar3 = bestParams.par3;
    elapsedtime = toc; 
 
end
