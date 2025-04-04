function [X_data, DX,t_ae, x_tn] = E_Data(model,g,tau,phi,par,tspan,t_r,n_a,N,dist)
%% Find solution
sol = dde23(@(t,y,Z) g(t,y,Z,par), tau, phi, tspan);

%% Select data points based on user input
if dist == 1
    t_ae = linspace(tspan(1), tspan(1) + ((tspan(2) - tspan(1))*t_r), N);
else
    t_ae = tspan(1) + ((tspan(2) - tspan(1))*t_r) * rand(1, N);
end

%%  Calculation of current and delayed state terms
numStates = size(sol.y, 1);
x_t = zeros(numStates, length(t_ae));
x_tn = zeros(numStates, length(t_ae));

for i = 1:numStates
    x_t(i, :) = deval(sol, t_ae, i);
    x_tn(i, :) = x_t(i, :) + n_a * randn(size(x_t(i, :)));
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
                x_d(j, k, i) = interp1(t_ae, x_tn(j, :), t_d(i), 'linear', 'extrap');
            end
        end
    end
end
%% derivative
DX = zeros(size(x_tn, 2), size(x_tn, 1));

for j = 1:size(x_tn, 2)
    Z = reshape(x_d(:, :, j), size(x_d, 1), size(x_d, 2));
    DX(j, :) = g(0, x_tn(:, j), Z, par).';
end

x_dn = reshape(x_d, numDelays*numStates, length(t_ae));

%% Construct non linearities for library
% if strcmp(model, 'MG')
%     X_data = [x_tn', x_dn', 1./ (1 + x_dn.^par(3))'];
% else
%       X_data = [x_tn', x_dn'];
% end
switch model
    case 'MG'
        X_data = [x_tn', x_dn', 1./ (1 + x_dn.^par(3))'];
    otherwise
        X_data = [x_tn', x_dn'];
end



