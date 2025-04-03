function [X_data,DX,Theta_all,bestTau,evalCount,elapsedtime,bestpar3] = BO(model,tau,par,DX_true, folderName) 


tic;   
Xi = [];
evalCount =0;
   switch model
        case {'Rossler', 'DLE','TDN'}
            bestTau=[];
        case 'MG'
            bestTau=[];
            bestPar3=[];
        otherwise
            bestTau=[];
    end


    switch model
        case 'Rossler'
            tauRange = optimizableVariable('tau', [1, 2.5], 'Transform','none');
        case {'DLE','TDN'}
            tauRange = optimizableVariable('tau', [0.1, 1.5], 'Transform','none');
        case 'MG'
            tauRange = optimizableVariable('tau', [0.1, 2], 'Transform','none');
            par3Range = optimizableVariable('par3', [0.1, 10], 'Transform','none');
        otherwise
            tauRange = optimizableVariable('tau', [0.1, 1.5], 'Transform','none');
    end

    

    function objective = optimizeParameters(params)
        switch model
            case {'Rossler', 'DLE','TDN'}
                [g,phi,par,tspan,X_data,DX,Xi,Theta_all,x_sindy,x_sol] = Popt(model, [tau(1),params.tau],par);
            case 'MG'
               [g,phi,par,tspan,X_data,DX,Xi,Theta_all,x_sindy,x_sol] = Popt(model, params.tau, [par(1),par(2),params.par3]);
            otherwise
                [g,phi,par,tspan,X_data,DX,Xi,Theta_all,x_sindy,x_sol] = Popt(model, params.tau,par);
        end
        err = norm(DX_true - Theta_all * Xi, 2)+norm(x_sol-x_sindy,2);
        objective = err; % Use norm error as the objective to minimize
       evalCount = evalCount + 1;
       
    end
    
    % Run Bayesian Optimization
    switch model
            case {'Rossler', 'DLE','TDN'}
                % results = bayesopt(@optimizeParameters, [tau1Range, tau2Range,lambdaRange, alphaRange], 'Verbose', 0, ...
                %        'MaxObjectiveEvaluations', 300);
                results = bayesopt(@optimizeParameters, tauRange, ...
        'IsObjectiveDeterministic', false, ...
        'AcquisitionFunctionName', 'expected-improvement-plus', ...
        'Verbose', 1, ...
        'ExplorationRatio', 0.6, ...
        'MaxObjectiveEvaluations', 300);%'TolFun',1e-6
            case 'MG'
                % results = bayesopt(@optimizeParameters, [tau1Range, tau2Range,lambdaRange, alphaRange], 'Verbose', 0, ...
                %        'MaxObjectiveEvaluations', 300);
                results = bayesopt(@optimizeParameters, [tauRange, par3Range], ...
        'IsObjectiveDeterministic', false, ...
        'AcquisitionFunctionName', 'expected-improvement-plus', ...
        'Verbose', 2, ...
        'ExplorationRatio', 0.6, ...
        'MaxObjectiveEvaluations', 300);%'TolFun',1e-6
            otherwise
%                 results = bayesopt(@optimizeParameters, [tauRange,lambdaRange, alphaRange], 'Verbose', 0, ...
%                        'MaxObjectiveEvaluations', 200);
                results = bayesopt(@optimizeParameters, tauRange, ...
        'IsObjectiveDeterministic', false, ...
        'AcquisitionFunctionName', 'expected-improvement-plus', ...
        'Verbose', 1, ...
        'ExplorationRatio', 0.6, ...
        'MaxObjectiveEvaluations', 300);
   
        end
    disp(results);
   disp(fieldnames(results));


    % Optimal lambda and alpha found by Bayesian Optimization
    switch model
        case {'Rossler', 'DLE','TDN'}
%             bestTau1 = results.XAtMinObjective.tau1;
%             bestTau2 = results.XAtMinObjective.tau2;
            bestTau = table2array(results.XAtMinEstimatedObjective(1, 'tau'));
        case 'MG'
            bestTau = table2array(results.XAtMinEstimatedObjective(1, 'tau'));
            bestPar3 = table2array(results.XAtMinEstimatedObjective(1, 'par3'));
        otherwise
%             bestTau = results.XAtMinObjective.tau;
            bestTau = table2array(results.XAtMinEstimatedObjective(1, 'tau'));
    end 
   switch model
        case {'Rossler', 'DLE','TDN'}
            bestTau=[tau(1),bestTau];
            case 'MG'
            bestTau=bestTau;
            bestPar=[par(1),par(2),bestPar3];
        otherwise
            bestTau = bestTau;
    end 
    
elapsedtime=toc;
    dataFileName = fullfile(folderName, 'BO_Optimization_Data.mat');
    save(dataFileName);
end
