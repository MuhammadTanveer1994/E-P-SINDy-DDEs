function [error_metrics] = Error(Theta_all, Theta_train, Theta_valid, Xi, DX, DX_train, DX_valid, T, t_ae, t_train, t_valid, solA, solSINDy)
 
    %% Complete Dataset Error Calculations
    t_all = linspace(0, T, 100);
    
    o_sol = deval(solA, t_all);        
    i_sindy = deval(solSINDy, t_all);    
    
    x_sol = deval(solA, t_ae);           
    x_sindy = deval(solSINDy, t_ae);
    
    a = Theta_all * Xi;
    
    rmse_complete_deriv = sqrt(mean((DX - a).^2, 'all'));
    rmse_complete_sample  = sqrt(mean((x_sol - x_sindy).^2, 'all'));
    rmse_complete_global  = sqrt(mean((o_sol - i_sindy).^2, 'all'));
    
    rmse_complete_combined = rmse_complete_deriv + rmse_complete_sample;
    
  
    %% Training Dataset Error Calculations
    a_train = Theta_train * Xi;
    
    x_sol_train = deval(solA, t_train);
    x_sindy_train = deval(solSINDy, t_train);
    
    rmse_train_deriv = sqrt(mean((DX_train - a_train).^2, 'all'));
    rmse_train_solution = sqrt(mean((x_sol_train - x_sindy_train).^2, 'all'));
    
    %% Validation Dataset Error Calculations
    a_valid = Theta_valid * Xi;
    
    x_sol_valid = deval(solA, t_valid);
    x_sindy_valid = deval(solSINDy, t_valid);
    
    rmse_valid_deriv = sqrt(mean((DX_valid - a_valid).^2, 'all'));
    rmse_valid_solution = sqrt(mean((x_sol_valid - x_sindy_valid).^2, 'all'));
    
    %% Organize Error Metrics into a Structure
    error_metrics.complete.rmse.derivatives      = rmse_complete_deriv;
    error_metrics.complete.rmse.solution_sample    = rmse_complete_sample;
    error_metrics.complete.rmse.solution_global    = rmse_complete_global;
    error_metrics.complete.rmse.combined           = rmse_complete_combined;

    error_metrics.training.rmse.derivatives      = rmse_train_deriv;
    error_metrics.training.rmse.solution         = rmse_train_solution;
    
    error_metrics.validation.rmse.derivatives    = rmse_valid_deriv;
    error_metrics.validation.rmse.solution       = rmse_valid_solution;
    
    %% Display the Error Metrics
    fprintf('\n===== Complete Dataset Error Metrics =====\n');
    fprintf('RMSE (Derivatives): %e\n', rmse_complete_deriv);
    fprintf('RMSE (Solution at Sample Points): %e\n', rmse_complete_sample);
    fprintf('RMSE (Global Solution): %e\n', rmse_complete_global);
    fprintf('RMSE (Combined Error): %e\n', rmse_complete_combined);

    fprintf('\n===== Training Dataset Error Metrics =====\n');
    fprintf('RMSE (Derivatives): %e\n', rmse_train_deriv);
    fprintf('RMSE (Solution): %e\n', rmse_train_solution);
    
    fprintf('\n===== Validation Dataset Error Metrics =====\n');
    fprintf('RMSE (Derivatives): %e\n', rmse_valid_deriv);
    fprintf('RMSE (Solution): %e\n', rmse_valid_solution);
    
end
