function [g,tau,phi,par,n] = P_Model(model)

 switch model
    case {'Rossler1','Rossler2'}
        phi =@(x) kron(ones(length(x),1),[1.5; 0.4; 0.9]);
        n = 3; % Rossler system

    case 'tau_3'
        phi =@(x) kron(ones(length(x),1),[0.5; -0.5]);
        n = 2; % This model has 2 state variables and three delays

    otherwise
        phi = @(x) cos(x);
        n = 1;
end
switch model
    case 'DL'  % Delay logistic equation
        g = @(t,y,Z,par) par(1) * y * (1 - Z/par(2));
        tau = 1;
        par = [1.8,10];

    case 'MG'   % Mackey glass equation
        g = @(t,y,Z,par) par(1) * Z / (1 + Z^par(3)) - par(2) * y ;
        tau = 1;
        par = [4, 2, 9.6];

    case 'Rossler1' % Rossler system with one delay
        g = @(t,y,Z,par) [-y(3) - y(2)+ par(1) * Z(1,1)+ par(2) * Z(1,2);...
            y(1) + par(3) * y(2); par(4) + y(3) * (y(1) - par(5))];
        tau = [0.0001, 2];
        par = [0, 1, 0.2, 0.2, 1.2];

    case 'Rossler2' % Rossler system with two delays
        g = @(t,y,Z,par) [-y(3) - y(2)+ par(1) * Z(1,1)+ par(2) * Z(1,2);...
            y(1) + par(3) * y(2); par(4) + y(3) * (y(1) - par(5))];
        tau = [1, 2];
        par = [0.2, 1, 0.2, 0.2, 1.2];

    case 'tau_3' % two neuron system
        g = @(t,y,Z,par) [-par(1)*y(1) + par(2)*tanh(Z(1,1))+  par(3)*tanh(Z(2,3));...
            -par(1)*y(2) + + par(2)*tanh(Z(2,1))+ par(4)*tanh(Z(1,2))];
        tau = [1.5, 2, 2];
        par = [0.5, -1, 1,2];
    otherwise
        error('Invalid model type.');
end
