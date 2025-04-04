function Xi = sparsifyDynamics(Theta_all, DX, lambda, alpha, n, method, k)
    Xi = [];
    
    switch method
        case 1
            % Sequentially thresholded least squares
            Xi = Theta_all \ DX; 
            for iter = 1:k 
                smallinds = (abs(Xi) < lambda);   
                Xi(smallinds) = 0;                
                for ind = 1:n                     
                    biginds = ~smallinds(:, ind);
                    Xi(biginds, ind) = Theta_all(:, biginds) \ DX(:, ind); 
                end
            end
        case 2
            % LASSO Regression
            Xi = zeros(size(Theta_all,2), n);
            for i = 1:n
                [B, FitInfo] = lasso(Theta_all, DX(:, i), 'Lambda', lambda, 'Alpha', alpha);
                Xi(:, i) = B;
            end
        case 3
               % Elastic Net Regression
                Xi = zeros(size(Theta_all, 2), n);
                for i = 1:n
                    [B, ~] = lasso(Theta_all, DX(:, i), 'Lambda', lambda, 'Alpha', alpha);
                    Xi(:, i) = B;
                end
    end
 
end

