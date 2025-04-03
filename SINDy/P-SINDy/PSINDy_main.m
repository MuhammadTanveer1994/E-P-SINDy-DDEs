% PSINDy_main.m - main script for running PSINDy analysis
clear all; close all; clc;

%% User inputs - Modify these values as needed
model = 'DL';        % Choose from: 'DL', 'MG', 'Rossler1', 'Rossler2', 'tau_3'
optimizer = 'none';    % Choose from: 'none', 'PSO' 'BO', 'BF'

%% PSINDy Parameters
polyorder = 2;       
usesine = 0;        
training_ratio = 0.8;
noise_amplitude = 0;
final_time = 32;
num_samples = 500;
M = 10;              % number of collocating polynomials

%% Sparse Regression Parameters
regression_method = 1;  % 1: Sequential Threshold, 2: LASSO
lambda = 1e-10;
alpha = 1e-10;
k = 20;

%% Create results folder with timestamp
folderName = sprintf('SimulationResults-one-delay-PSINDy-%s-%s-M%d', model, optimizer, M);

if ~exist(folderName, 'dir')
    mkdir(folderName);
end

% Save input parameters
inputParams = struct();
inputParams.model = model;
inputParams.optimizer = optimizer;
inputParams.polyorder = polyorder;
inputParams.usesine = usesine;
inputParams.training_ratio = training_ratio;
inputParams.noise_amplitude = noise_amplitude;
inputParams.final_time = final_time;
inputParams.num_samples = num_samples;
inputParams.M = M;
inputParams.regression_method = regression_method;
inputParams.lambda = lambda;
inputParams.alpha = alpha;
inputParams.k = k;

save(fullfile(folderName, 'input_parameters.mat'), '-struct', 'inputParams');

%% Start timing
tic;

%% Get model details
[g, tau, phi, par, n] = P_Model(model);
tspan = [0 final_time];

%% Generate data
[X_data, DX_true, u0, D, t_ae, x_tn] = P_Data(model, g, tau, phi, par, tspan, ...
    M, training_ratio, noise_amplitude, num_samples, 1); % Using uniform distribution

%% Split data into training and validation sets
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
            bestpar3 = par(3);
        end

    case 'BO'
                 if strcmp(model, 'MG')
        [X_data, DX, Theta_all, bestTau, evalCount, elapsedTime,bestpar3] = BO(model, tau, par, DX_true, folderName);
                 else
        [X_data, DX, Theta_all, bestTau, evalCount, elapsedTime] = BO(model, tau, par, DX_true, folderName);

                 end
        optimization_results.method = 'BO';
        optimization_results.evalCount = evalCount;
        optimization_results.elapsedTime = elapsedTime;
        
    case 'PSO'
        if strcmp(model, 'MG')
        [bestTau, evalCount, elapsedTime, Theta_all, DX, X_data,bestpar3] = PSO(model, tau, par, DX_true, folderName);
        else
        [bestTau, evalCount, elapsedTime, Theta_all, DX, X_data] = PSO(model, tau, par, DX_true, folderName);
        end
        optimization_results.method = 'PSO';
        optimization_results.evalCount = evalCount;
        optimization_results.elapsedTime = elapsedTime;
      
    case 'BF'
         if strcmp(model, 'MG')
        [bestTau, evalCount, elapsedTime, Theta_all, DX,bestpar3] = BF(model, tau, DX_true, folderName);
         else    
        [bestTau, evalCount, elapsedTime, Theta_all, DX] = BF(model, tau, DX_true, folderName);
         end
        optimization_results.method = 'BF';
        optimization_results.evalCount = evalCount;
        optimization_results.elapsedTime = elapsedTime;
end

save(fullfile(folderName, 'optimization_results.mat'), '-struct', 'optimization_results');

if ~strcmp(optimizer, 'none') && ~strcmp(optimizer, 'BF')
    train_size = round(training_ratio * length(t_ae));
    t_train = t_ae(1:train_size);
    t_valid = t_ae(train_size+1:end);
    X_train = X_data(1:train_size, :);
    X_valid = X_data(train_size+1:end, :);
    DX_train = DX(1:train_size, :);
    DX_valid = DX(train_size+1:end, :);
    
    Theta_train = poolData(X_train, polyorder, usesine);
    Theta_valid = poolData(X_valid, polyorder, usesine);
end

%% Compute sparse regression
Xi = sparsifyDynamics(Theta_all, DX, lambda, alpha, n, regression_method, k);
disp(Xi)

m = size(X_data, 2);
library(m, polyorder);

%% Validation
ode_options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
dde_options = ddeset('RelTol', 1e-5, 'AbsTol', 1e-5, 'MaxStep', 1e-2);

soltrue = dde23(@(t,y,Z) g(t,y,Z,par), tau, phi, tspan, dde_options);

bestTau_new = max(abs(bestTau));
[D_new, theta_new] = difmat(-bestTau_new, 0, M);
X_new = phi(theta_new);
u0_new = X_new;

solA = ode45(@(t,u) d_ODE(t, u, par, tau, M, model, D), tspan, u0, ode_options);
 if strcmp(model, 'MG')
solSINDy = ode45(@(t,u) gPSINDy(t, u, Xi, polyorder, usesine, model, D_new,bestpar3), tspan, u0_new, ode_options);
 else
     solSINDy = ode45(@(t,u) gPSINDy(t, u, Xi, polyorder, usesine, model, D_new), tspan, u0_new, ode_options);
 end
%% Plot results
plotResult( solA, solSINDy, training_ratio, t_ae, tspan, x_tn, ...
    final_time, num_samples, M, noise_amplitude, model, 1, n, folderName);

%% Calculate errors using the enhanced Error function
[error_metrics] = Error(Theta_all, Theta_train, Theta_valid, Xi, DX, DX_train, DX_valid, final_time, t_ae, t_train, t_valid, solA, solSINDy);

%% Save simulation results
simulation_results = struct();
simulation_results.model_details.g = g;
simulation_results.model_details.tau = tau;
simulation_results.model_details.bestTau = bestTau;
simulation_results.model_details.phi = phi;
simulation_results.model_details.par = par;
simulation_results.model_details.tspan = tspan;
simulation_results.model_details.M = M;

simulation_results.data.X_data = X_data;
simulation_results.data.DX_true = DX_true;
simulation_results.data.DX = DX;
simulation_results.data.t_ae = t_ae;
simulation_results.data.x_tn = x_tn;
simulation_results.data.u0 = u0;
simulation_results.data.D = D;
simulation_results.data.D_new = D_new;
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

simulation_results.solutions.solA = solA;
simulation_results.solutions.solSINDy = solSINDy;

simulation_results.errors = error_metrics;

simulation_results.timing.total_time = toc;

% Save complete results
save(fullfile(folderName, 'simulation_results.mat'), '-struct', 'simulation_results');

% Save a summary text file
fid = fopen(fullfile(folderName, 'simulation_summary.txt'), 'w');
fprintf(fid, 'PSINDy Simulation Summary\n');
fprintf(fid, '========================\n\n');

fprintf(fid, 'Model: %s\n', model);
fprintf(fid, 'Optimizer: %s\n', optimizer);
fprintf(fid, 'Collocation Points (M): %d\n\n', M);

fprintf(fid, 'Error Metrics:\n');
fprintf(fid, 'Complete Data:\n');
fprintf(fid, '  RMSE (Derivatives): %e\n', error_metrics.complete.rmse.derivatives);
fprintf(fid, '  RMSE (Solution at Sample Points): %e\n', error_metrics.complete.rmse.solution_sample);

fprintf(fid, '\nTraining Data:\n');
fprintf(fid, '  RMSE (Derivatives): %e\n', error_metrics.training.rmse.derivatives);
fprintf(fid, '  RMSE (Solution): %e\n', error_metrics.training.rmse.solution);

fprintf(fid, '\nValidation Data:\n');
fprintf(fid, '  RMSE (Derivatives): %e\n', error_metrics.validation.rmse.derivatives);
fprintf(fid, '  RMSE (Solution): %e\n', error_metrics.validation.rmse.solution);


fprintf(fid, '\nTotal Execution Time: %.2f seconds\n', toc);
fclose(fid);

%% Display results summary
fprintf('\nPSINDy Simulation completed successfully!\n');
fprintf('Results saved in folder: %s\n', folderName);
fprintf('Files saved:\n');
fprintf('1. input_parameters.mat - All input parameters\n');
fprintf('2. optimization_results.mat - Optimization details\n');
fprintf('3. simulation_results.mat - Complete simulation results\n');
fprintf('4. simulation_summary.txt - Text summary of results\n');
fprintf('Total execution time: %.2f seconds\n', toc);