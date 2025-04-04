function [X_data,DX,u0,D,t_ae, x_tn] = P_Data(model,g,tau,phi,par,tspan,M,t_r,n_a,N,dist)
     %% Find solution
    soltrue = dde23(@(t,y,Z) g(t,y,Z,par), tau, phi, tspan);
  
    %% Select data points based on user input
    if dist == 1
        t_ae = linspace(tspan(1), tspan(1) + ((tspan(2) - tspan(1))*t_r), N);
    else
        t_ae = tspan(1) + ((tspan(2) - tspan(1))*t_r) * rand(1, N);
    end
   
   %%  Calculation of current and delayed state terms 
numStates = size(soltrue.y, 1); 
x_to = zeros(numStates, length(t_ae));
x_ton = zeros(numStates, length(t_ae));

for i = 1:numStates
    x_to(i, :) = deval(soltrue, t_ae, i);
    x_ton(i, :) = x_to(i, :) + n_a * randn(size(x_to(i, :)));
end

numDelays = size(tau,2);

x_d = zeros(numStates, numDelays, length(t_ae)); 
for k = 1:numDelays
    t_d = t_ae - tau(k);
    for i = 1:length(t_d)
        if t_d(i) <= 0
            x_d(:, k, i) = phi(t_d(i)); 
        else
            for j = 1:numStates
                x_d(j, k, i) = interp1(t_ae, x_ton(j, :), t_d(i), 'linear', 'extrap');
            end
        end
    end
end
%% derivative
DX_o = zeros(size(x_ton, 2), size(x_ton, 1)); 

for j = 1:size(x_ton, 2)  
    Z = reshape(x_d(:, :, j), size(x_d, 1), size(x_d, 2)); 
    DX_o(j, :) = g(0, x_ton(:, j), Z, par).'; 
end

    %% Calculations for library construction
    tau_max=max(abs(tau));
    [D, theta] = difmat(-tau_max, 0,M);
    X = phi(theta);
    u0=X;
    
    options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
    sol = ode45(@(t,u) d_ODE(t, u, par, tau, M, model,D), tspan, u0,options);

    %% Find current and delayed state terms
    x_t = deval(sol, t_ae);
    x_tn = x_t + n_a * randn(size(x_t));

%% Compute the derivative at each time point
for i = 1:length(x_tn)
    t = t_ae(i);
    u = deval(sol, t); 
    DX(:,i)=   d_ODE(0, u, par, tau, M, model,D);
end

DX=DX';

%% data
    switch model
        case 'MG'
            X_data = [x_tn', 1./ (1 + x_tn.^par(3))'];
        otherwise
             X_data = x_tn';
    end

end



