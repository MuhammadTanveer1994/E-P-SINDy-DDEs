function [g,tau,phi,par,tspan,sol,X_data, DX_true,Xi,Theta_all,solSINDy] = ESINDy(model)

folderName = sprintf('SimulationResults-SQTL-PSO1-100-%s', num2str(model)); % Adjust model to your naming scheme if needed
if ~exist(folderName, 'dir')
    mkdir(folderName); % Create the folder if it doesn't exist
end
  tic;
%% DDE Model
    [g,tau,phi,par,n] = E_Model(model);
    
    %% SINDy Parameters
    polyorder=2; 
    usesine=1;
    t_r=0.6; % training ratio
    n_a = 0; % noise amplitude
    T=30; tspan=[0 T];  % T final time 
    N=100;% 'number of samples'

    %% Compute Derivative and Construct the Data Matrix
    dist = 1;   % '(1) Uniform', '(2) Random'
    [X_data, DX_true, t_ae, x_tn] = E_Data(model, g, tau, phi, par,tspan, t_r, n_a, N, dist);
     
    %% Build Library with Delay
    Theta_all = poolData(X_data, polyorder,usesine);
% DX=DX_true;
    %% Sparse Regression
    method = 1; % '1. Sequentially Thresholded Least Square Regression', '2. LASSO Regression',, '3. Elastic Net Regression'
    lambda = 1e-15;  alpha = 1e-15; k=20;
     % [X_data,DX,Theta_all,bestTau,evalCount,elapsedtime] = BO(model,tau,par,DX_true);
    % [Xi, bestTau, evalCount,elapsedtime,Theta_all,DX] = BF(model,tau,DX_true,folderName);
    % [bestTau,evalCount,elapsedtime,Theta_all,DX] = PSO(model,tau,par,DX_true);
    % [bestTau,evalCount,elapsedtime,Theta_all,DX,alpha,lambda,k] = PSO(model,tau,par,DX_true);

    %%
%[bestTau, evalCount, elapsedtime, Theta_all, DX] = PSO2(model, tau, par, DX_true, folderName);



  [bestTau, evalCount, elapsedtime, Theta_all, DX] = PSO1(model, tau,par, DX_true, folderName); % only optimize tau  
    
    % [bestTau, evalCount, elapsedtime, Theta_all, DX, alpha, lambda, k] = PSO(model, tau, par, DX_true, folderName); %optimize hyperparameters and tau
    % [Xi, bestLambda, bestK, bestAlpha] = sparsifyDynamics_opt(Theta_all, DX, n, method, folderName);%BF
    % [Xi, bestLambda, bestAlpha, bestK] = sparsifyDynamics_opt_PSO(Theta_all, DX, n, method)
    % Xi = sparsifyDynamics(Theta_all, DX, lambda, alpha, n,method);

%     % For PSO optimization (default)
% [bestTau, evalCount, elapsedtime, Theta_all, DX] = OptimizeWithPSOandBSA(model, tau, par, DX_true, folderName);

% For BSA optimization
%[bestTau, evalCount, elapsedtime, Theta_all, DX] = OptimizeWithPSOandBSA(model, tau, par, DX_true, folderName, 'BSA');
    Xi = sparsifyDynamics(Theta_all, DX, lambda, alpha, n, method, k);
    m=size(X_data,2);
    library(m,polyorder)
    % bestTau=tau;
    %% Validation
    options = ddeset('RelTol', 1e-5, 'AbsTol', 1e-5, 'MaxStep', 1e-2);
    % options = ddeset('RelTol', 1e-7, 'AbsTol', 1e-7);
    sol = dde23(@(t,y, Z) g(t,y, Z, par), tau, phi, tspan,options);
    solSINDy = dde23(@(t,y, Z) gSINDy(t, y, Z, Xi, polyorder,usesine, model), bestTau, phi, tspan, options);
    
    %% Plotting
    
    plotResults(sol, solSINDy,t_r,t_ae, x_tn,T,N,n_a,model,dist,folderName);
   
    %%  Errors calculation
    [err_1,err_2a,err_2b,err_3]=Error(Theta_all,Xi,DX,T,t_ae,sol,solSINDy);
       elapsedtime = toc;

filename = fullfile(folderName, sprintf('ESINDySimulationData-%s.mat', num2str(model)));

% filename = sprintf('ESINDySimulationData-%s.mat', num2str(model));
save(filename);
end
