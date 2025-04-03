function [g, phi, par, tspan, X_data, DX, Xi, Theta_all, x_sindy, x_sol] = Eopt(model, tauVal, parVal)

[g,~,phi,par,n] = E_Model(model);
if strcmp(model, 'MG')
        par(3) = parVal;      
        bestPar3 = par(3);
    end

polyorder=2;
usesine=0;
t_r=0.6;
n_a = 0;
T=30; tspan=[0 T];
N=100;
lambda = 1e-10;  alpha = 1e-10; k=20;
dist = 1;   % '(1) Uniform', '(2) Random'
[X_data, DX, t_ae, x_tn] = E_Data(model, g, tauVal, phi, par, tspan, t_r, n_a, N, dist);

Theta_all = poolData(X_data, polyorder,usesine);
tau=tauVal;
%% Sparse Regression
method = 1; % '1. Sequentially Thresholded Least Square Regression', '2. LASSO Regression'

Xi = sparsifyDynamics(Theta_all, DX, lambda, alpha, n, method, k);

options = ddeset('RelTol', 1e-5, 'AbsTol', 1e-5, 'MaxStep', 1e-2);
sol = dde23(@(t,y, Z) g(t,y, Z, par), tau, phi, tspan,options);
 if strcmp(model, 'MG')
solSINDy = dde23(@(t,y, Z) gSINDy(t, y, Z, Xi, polyorder,usesine, model,bestPar3), tau, phi, tspan, options);
 else
solSINDy = dde23(@(t,y, Z) gSINDy(t, y, Z, Xi, polyorder,usesine, model), tau, phi, tspan, options);
 end
t_all = linspace(0,T,N);
x_sindy=deval(solSINDy,t_ae);
x_sol = deval(sol, t_ae);
end
