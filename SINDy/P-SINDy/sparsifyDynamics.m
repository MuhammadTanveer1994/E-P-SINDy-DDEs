function Xi = sparsifyDynamics(Theta_all, DX, lambda, alpha, n, method, k)
    Xi = [];

    if method == 1
        % Sequentially thresholded least squares
        Xi = Theta_all \ DX; 
        for k = 1:k
            smallinds = (abs(Xi) < lambda);   
            Xi(smallinds) = 0;                
            for ind = 1:n                    
                biginds = ~smallinds(:, ind);
                Xi(biginds, ind) = Theta_all(:, biginds) \ DX(:, ind); 

            end
        end
    elseif method == 2
        % LASSO Regression
        Xi = zeros(size(Theta_all,2), n);
        for i = 1:n
            [B, FitInfo] = lasso(Theta_all, DX(:, i), 'Lambda', lambda, 'Alpha', alpha);
            Xi(:, i) = B;
        end
    end
 
end
