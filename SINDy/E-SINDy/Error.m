function [error_metrics] = Error(Theta_all, Theta_train, Theta_valid, Xi, DX, DX_train, DX_valid, T, t_ae, t_train, t_valid, sol, solSINDy)
    error_metrics = struct();
    
    %% Complete dataset errors
    a = Theta_all * Xi;
    
    t_all = linspace(0, T, 100);
    o_sol = deval(sol, t_all);        
    i_sindy = deval(solSINDy, t_all); 
    x_sindy = deval(solSINDy, t_ae);  
    x_sol = deval(sol, t_ae);         
    
    
    rmse_1 = sqrt(mean((DX - a).^2, 'all'));         
    rmse_2a = sqrt(mean((x_sol - x_sindy).^2, 'all')); 
    rmse_2b = sqrt(mean((o_sol - i_sindy).^2, 'all')); 
    
    error_metrics.complete.rmse.derivatives = rmse_1;
    error_metrics.complete.rmse.solution_sample = rmse_2a;
    error_metrics.complete.rmse.solution_full = rmse_2b;
    
    %% Training dataset errors
    a_train = Theta_train * Xi;
    
    x_sindy_train = deval(solSINDy, t_train);
    x_sol_train = deval(sol, t_train);
    
    rmse_train_1 = sqrt(mean((DX_train - a_train).^2, 'all'));  
    rmse_train_2 = sqrt(mean((x_sol_train - x_sindy_train).^2, 'all')); 
    
    error_metrics.training.rmse.derivatives = rmse_train_1;
    error_metrics.training.rmse.solution = rmse_train_2;
    
    %% Validation dataset errors
    a_valid = Theta_valid * Xi;
    
    x_sindy_valid = deval(solSINDy, t_valid);
    x_sol_valid = deval(sol, t_valid);
        
    rmse_valid_1 = sqrt(mean((DX_valid - a_valid).^2, 'all'));  
    rmse_valid_2 = sqrt(mean((x_sol_valid - x_sindy_valid).^2, 'all')); 
    
    error_metrics.validation.rmse.derivatives = rmse_valid_1;
    error_metrics.validation.rmse.solution = rmse_valid_2;
        
    fprintf('\n===== Error Metrics =====\n');
    fprintf('Complete Dataset:\n');
    fprintf('  RMSE (Derivatives): %e\n', rmse_1);
    fprintf('  RMSE (Solution at Sample Points): %e\n', rmse_2a);
    fprintf('  RMSE (Full Solution): %e\n', rmse_2b);
    
    fprintf('\nTraining Dataset:\n');
    fprintf('  RMSE (Derivatives): %e\n', rmse_train_1);
    fprintf('  RMSE (Solution): %e\n', rmse_train_2);
    
    fprintf('\nValidation Dataset:\n');
    fprintf('  RMSE (Derivatives): %e\n', rmse_valid_1);
    fprintf('  RMSE (Solution): %e\n', rmse_valid_2);
    
end