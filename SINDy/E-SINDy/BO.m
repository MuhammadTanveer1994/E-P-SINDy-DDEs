function [X_data,DX,Theta_all,bestTau,evalCount,elapsedtime,bestPar3] = BO(model,tau,par,DX_true, folderName) 

tic;   
Xi = [];
evalCount =0;
    switch model
        case 'tau_3'
            bestTau=[];
        case {'Rossler1', 'Rossler2'}           
            bestTau=[];
        case 'MG'
            bestTau=[];
            bestPar3=[];
        otherwise
            bestTau=[];
    end
    % Define the range
    switch model
        case 'tau_3'
            tau1Range = optimizableVariable('tau1', [1, 1.8], 'Transform','none');
            tau2Range = optimizableVariable('tau2', [1.5, 2.5], 'Transform','none');
            tau3Range = optimizableVariable('tau3', [1.5, 2.5], 'Transform','none');
        case 'Rossler2'
            tau1Range = optimizableVariable('tau1', [0.1, 1], 'Transform','none');
            tau2Range = optimizableVariable('tau2', [1, 3], 'Transform','none');
            case 'Rossler1'
            tau1Range = optimizableVariable('tau1', [0.000001, 1], 'Transform','none');
            tau2Range = optimizableVariable('tau2', [1, 3], 'Transform','none');
        case 'MG'
            tauRange = optimizableVariable('tau', [0.1, 2], 'Transform','none');
            par3Range = optimizableVariable('par3', [0.1, 20], 'Transform','none');
        otherwise
            tauRange = optimizableVariable('tau', [0.1, 1.5], 'Transform','none');
    end

    

    function objective = optimizeParameters(params)
        switch model
             case 'tau_3'
                [g,phi,par,tspan,X_data, DX,Xi,Theta_all,x_sindy,x_sol] = Eopt(model, [params.tau1, params.tau2,params.tau3]);
            case {'Rossler1', 'Rossler2'}
                [g,phi,par,tspan,X_data, DX,Xi,Theta_all,x_sindy,x_sol] = Eopt(model, [params.tau1, params.tau2]);
            case 'MG'
               [g,phi,par,tspan,X_data, DX,Xi,Theta_all,x_sindy,x_sol] = Eopt(model, params.tau, [params.par3]);
            otherwise
                [g,phi,par,tspan,X_data, DX,Xi,Theta_all,x_sindy,x_sol]  = Eopt(model, params.tau);
        end
        err = norm(DX_true - Theta_all * Xi, 2)+norm(x_sol-x_sindy,2);
        objective = err; 
       evalCount = evalCount + 1;
       
    end
    
    switch model
                case 'tau_3'
                results = bayesopt(@optimizeParameters, [tau1Range, tau2Range,tau3Range], ...
        'IsObjectiveDeterministic', false, ...
        'AcquisitionFunctionName', 'expected-improvement-plus', ...
        'Verbose', 1, ...
        'ExplorationRatio', 0.6, ...
        'MaxObjectiveEvaluations', 300);%'TolFun',1e-6
            case {'Rossler1', 'Rossler2'}
                results = bayesopt(@optimizeParameters, [tau1Range, tau2Range], ...
        'IsObjectiveDeterministic', false, ...
        'AcquisitionFunctionName', 'expected-improvement-plus', ...
        'Verbose', 1, ...
        'ExplorationRatio', 0.6, ...
        'MaxObjectiveEvaluations', 300);%'TolFun',1e-6
            case 'MG'
                results = bayesopt(@optimizeParameters, [tauRange, par3Range], ...
        'IsObjectiveDeterministic', false, ...
        'AcquisitionFunctionName', 'expected-improvement-plus', ...
        'Verbose', 1, ...
        'ExplorationRatio', 0.6, ...
        'MaxObjectiveEvaluations', 300);%'TolFun',1e-6
            otherwise
                results = bayesopt(@optimizeParameters, tauRange, ...
        'IsObjectiveDeterministic', false, ...
        'AcquisitionFunctionName', 'expected-improvement-plus', ...
        'Verbose', 0, ...
        'ExplorationRatio', 0.6, ...
        'MaxObjectiveEvaluations', 300,'OutputFcn',@terminationFunc);
   
        end
    disp(results);
   disp(fieldnames(results));


    switch model
         case 'tau_3'
            bestTau1 = table2array(results.XAtMinEstimatedObjective(1, 'tau1'));
            bestTau2 = table2array(results.XAtMinEstimatedObjective(1, 'tau2'));
             bestTau3 = table2array(results.XAtMinEstimatedObjective(1, 'tau3'));
        case {'Rossler1', 'Rossler2'}
            bestTau1 = table2array(results.XAtMinEstimatedObjective(1, 'tau1'));
            bestTau2 = table2array(results.XAtMinEstimatedObjective(1, 'tau2'));
        case 'MG'
            bestTau = table2array(results.XAtMinEstimatedObjective(1, 'tau'));
            bestPar3 = table2array(results.XAtMinEstimatedObjective(1, 'par3'));
        otherwise
            bestTau = table2array(results.XAtMinEstimatedObjective(1, 'tau'));
    end 
   switch model
       case 'tau_3'
            bestTau=[bestTau1,bestTau2,bestTau3];
        case {'Rossler1', 'Rossler2'}
            bestTau=[bestTau1,bestTau2];
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

