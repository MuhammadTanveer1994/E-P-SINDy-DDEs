% ESINDy_main.m - main script for running ESINDy analysis
clear all; close all; clc;

%% User inputs - Modify these values as needed
model = 'Rossler2';        % Choose from: 'DL', 'MG','Rossler1','Rossler2', 'tau_3'
optimizer = 'PSO';  % Choose from: 'none', 'PSO','BO', 'BF'

%% SINDy Parameters
polyorder = 2;
usesine = 0;
training_ratio = 0.6;
noise_amplitude = 0;
final_time = 30;
num_samples = 100;
distribution = 1;    % 1 for Uniform, 2 for Random

%% Sparse Regression Parameters
regression_method = 1;  % 1: Sequential Threshold, 2: LASSO
lambda = 1e-10;
alpha = 1e-10;
k = 20;

%% Create results folder
folderName = sprintf('SimulationResults-test-%s-%s-%d', model, optimizer, num_samples);
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

inputParams = struct();
inputParams.model = model;
inputParams.optimizer = optimizer;
inputParams.polyorder = polyorder;
inputParams.usesine = usesine;
inputParams.training_ratio = training_ratio;
inputParams.noise_amplitude = noise_amplitude;
inputParams.final_time = final_time;
inputParams.num_samples = num_samples;
inputParams.distribution = distribution;
inputParams.regression_method = regression_method;
inputParams.lambda = lambda;
inputParams.alpha = alpha;
inputParams.k = k;

save(fullfile(folderName, 'input_parameters.mat'), '-struct', 'inputParams');

%% Start timing
tic;

%% Get model details
[g, tau, phi, par, n] = E_Model(model);
tspan = [0 final_time];
%% Generate data
[X_data, DX_true, t_ae, x_tn] = E_Data(model, g, tau, phi, par, tspan,training_ratio, noise_amplitude, num_samples, distribution);

% Split data into training and validation sets
train_size = round(training_ratio * length(t_ae));
t_train = t_ae(1:train_size);
t_valid = t_ae(train_size+1:end);
X_train = X_data(1:train_size, :);
X_valid = X_data(train_size+1:end, :);
DX_train = DX_true(1:train_size, :);
DX_valid = DX_true(train_size+1:end, :);

%% Build library
Theta_all = poolData(X_data, polyorder, usesine);
Theta_train = poolData(X_train, polyorder, usesine);
Theta_valid = poolData(X_valid, polyorder, usesine);

%% Apply optimization if selected
optimization_results = struct();
switch optimizer
    case 'none'
        bestTau = tau;
        DX = DX_true;
        optimization_results.method = 'none';
        if strcmp(model, 'MG')
            optPar3 = par(3);
        end

    case 'PSO'
        if strcmp(model, 'MG')
        [bestTau, evalCount, elapsedTime, Theta_all, DX,optPar3] = PSO(model, tau, par, DX_true, folderName);
        optimization_results.optimized_par3 = optPar3;
        else
          [bestTau, evalCount, elapsedTime, Theta_all, DX] = PSO(model, tau, par, DX_true, folderName);

        end
        optimization_results.method = 'PSO';
        optimization_results.evalCount = evalCount;
        optimization_results.elapsedTime = elapsedTime;
       
        
        
    case 'BO'
        if strcmp(model, 'MG')
        [X_data, DX, Theta_all, bestTau, evalCount, elapsedTime,optPar3] = BO(model, tau, par, DX_true, folderName);
        optimization_results.optimized_par3 = optPar3;
        else
          [X_data, DX, Theta_all, bestTau, evalCount, elapsedTime] = BO(model, tau, par, DX_true, folderName);

        end
        optimization_results.method = 'BO';
        optimization_results.evalCount = evalCount;
        optimization_results.elapsedTime = elapsedTime;
       
        
    case 'BF'
    if strcmp(model, 'MG')
        [bestTau, evalCount, elapsedTime, Theta_all, DX, bestData, optPar3] = BF(model, tau, DX_true, par, folderName);
        optimization_results.optimized_par3 = optPar3;
    else
        [bestTau, evalCount, elapsedTime, Theta_all, DX] = BF(model, tau, DX_true, par, folderName);
    end            
        optimization_results.method = 'BF';
        optimization_results.evalCount = evalCount;
        optimization_results.elapsedTime = elapsedTime;

end

save(fullfile(folderName, 'optimization_results.mat'), '-struct', 'optimization_results');

%% Compute sparse regression
Xi = sparsifyDynamics(Theta_all, DX, lambda, alpha, n, regression_method, k);
disp(Xi)

%% library
m = size(X_data, 2);
library(m, polyorder);

%% Validation
options = ddeset('RelTol',1e-5,'AbsTol',1e-5,'MaxStep',1e-2);
sol = dde23(@(t,y,Z) g(t,y,Z,par), tau, phi, tspan, options);
if strcmp(model, 'MG')

    solSINDy = dde23(@(t,y,Z) gSINDy(t,y,Z,Xi,polyorder,usesine,model,optPar3), bestTau, phi, tspan, options);
else
    solSINDy = dde23(@(t,y,Z) gSINDy(t,y,Z,Xi,polyorder,usesine,model), bestTau, phi, tspan, options);
end

%% Plot results
plotResults(sol, solSINDy, training_ratio, t_ae, x_tn, final_time, ...
    num_samples, noise_amplitude, model, distribution, folderName);

%% Calculate errors 
[error_metrics] = Error(Theta_all, Theta_train, Theta_valid, Xi, DX_true, DX_train, DX_valid, final_time, t_ae, t_train, t_valid, sol, solSINDy);

%% Save all simulation results
simulation_results = struct();
simulation_results.model_details.g = g;
simulation_results.model_details.tau = tau;
simulation_results.model_details.bestTau = bestTau;
simulation_results.model_details.phi = phi;
simulation_results.model_details.par = par;
simulation_results.model_details.tspan = tspan;

simulation_results.data.X_data = X_data;
simulation_results.data.DX_true = DX_true;
simulation_results.data.DX = DX;
simulation_results.data.t_ae = t_ae;
simulation_results.data.x_tn = x_tn;
simulation_results.data.X_train = X_train;
simulation_results.data.X_valid = X_valid;
simulation_results.data.DX_train = DX_train;
simulation_results.data.DX_valid = DX_valid;
simulation_results.data.t_train = t_train;
simulation_results.data.t_valid = t_valid;

simulation_results.sindy.Xi = Xi;
simulation_results.sindy.Theta_all = Theta_all;
simulation_results.sindy.Theta_train = Theta_train;
simulation_results.sindy.Theta_valid = Theta_valid;

simulation_results.solutions.sol = sol;
simulation_results.solutions.solSINDy = solSINDy;

simulation_results.errors = error_metrics;

simulation_results.timing.total_time = toc;

save(fullfile(folderName, 'simulation_results.mat'), '-struct', 'simulation_results');

fid = fopen(fullfile(folderName, 'simulation_summary.txt'), 'w');
fprintf(fid, 'ESINDy Simulation Summary\n');
fprintf(fid, '========================\n\n');
fprintf(fid, 'Model: %s\n', model);
fprintf(fid, 'Optimizer: %s\n\n', optimizer);
fprintf(fid, 'Error Metrics:\n');
fprintf(fid, 'Complete Data:\n');
fprintf(fid, '  RMSE (Derivatives): %e\n', error_metrics.complete.rmse.derivatives);
fprintf(fid, '  RMSE (Solution at Sample Points): %e\n', error_metrics.complete.rmse.solution_sample);
fprintf(fid, '  RMSE (Full Solution): %e\n', error_metrics.complete.rmse.solution_full);
fprintf(fid, '\nTraining Data:\n');
fprintf(fid, '  RMSE (Derivatives): %e\n', error_metrics.training.rmse.derivatives);
fprintf(fid, '  RMSE (Solution): %e\n', error_metrics.training.rmse.solution);
fprintf(fid, '\nValidation Data:\n');
fprintf(fid, '  RMSE (Derivatives): %e\n', error_metrics.validation.rmse.derivatives);
fprintf(fid, '  RMSE (Solution): %e\n', error_metrics.validation.rmse.solution);
fprintf(fid, '\nTotal Execution Time: %.2f seconds\n', toc);
fclose(fid);

%% Display results summary
fprintf('\nSimulation completed successfully!\n');
fprintf('Results saved in folder: %s\n', folderName);
fprintf('Files saved:\n');
fprintf('1. input_parameters.mat - All input parameters\n');
fprintf('2. optimization_results.mat - Optimization details\n');
fprintf('3. simulation_results.mat - Complete simulation results\n');
fprintf('4. simulation_summary.txt - Text summary of results\n');
fprintf('Total execution time: %.2f seconds\n', toc);
