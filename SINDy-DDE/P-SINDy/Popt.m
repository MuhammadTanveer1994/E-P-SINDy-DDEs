function [g,phi,par,tspan,X_data,DX,Xi,Theta_all] = Popt(model,tau,par)

if strcmp(model, 'MG')
        [g,~,phi,par0,n] = P_Model(model);
        par0(3) = par(3);  
        par = par0;
    else
        [g,~,phi,par,n] = P_Model(model);
    end

    %% Parameters 
    polyorder = 2;   
    usesine = 1;
    t_r = 0.6; 
    n_a = 0; 
    T = 30; 
    tspan = [0 T]; 
    N = 100;
    M = 15;  

    dist = 1;  % '(1) Uniform', '(2) Random';
    [X_data, DX, u0, D, t_ae, x_tn] = P_Data(model, g, tau, phi, par, tspan, M, t_r, n_a, N, dist);

    Theta_all = poolData(X_data, polyorder, usesine);

    method = 1; % '1. STLS', '2. LASSO'
    lambda = 1e-10;  
    alpha = 1e-10;
    k = 20;

    Xi = sparsifyDynamics(Theta_all, DX, lambda, alpha, n, method, k);

    % options = odeset('RelTol',1e-4,'AbsTol',1e-6);
    % 
    % bestTau_new = max(abs(tau));
    % [D_new, theta_new] = difmat(-bestTau_new, 0, M);
    % X_new = phi(theta_new);
    % u0_new = X_new;
    % 
    % solA = ode15s(@(t,u) d_ODE(t, u, par, tau, M, model, D), tspan, u0, options);
    % 
    % if strcmp(model, 'MG')
    %    solSINDy = ode15s(@(t,u) gPSINDy(t, u, Xi, polyorder, usesine, model, D_new, par(3)), tspan, u0_new, options);
    % else
    %    solSINDy = ode15s(@(t,u) gPSINDy(t, u, Xi, polyorder, usesine, model, D_new), tspan, u0_new, options);
    % end
    % 
    % t_valid = solSINDy.x;  
    % t_ae = t_ae(t_ae >= t_valid(1) & t_ae <= t_valid(end));  
    % 
    % x_sindy = deval(solSINDy, t_ae);
    % x_sol = deval(solA, t_ae);
end
