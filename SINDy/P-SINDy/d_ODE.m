function dydt = d_ODE(t,u, par, tau, M, model,D)

switch model
    case {'Rossler1', 'Rossler2'}
        tau_max=max(abs(tau));
        [D,theta]=difmat(-tau_max,0,M);
        w=barywei(theta);
        d2=3;
        dMDM_DDE = kron(D(2:end,:), eye(d2));
        PM=@(u,j,x)barint(theta,w,u(j:d2:end)',x);

    case 'tau_3'
        tau_max=max(abs(tau));
        [D,theta]=difmat(-tau_max,0,M);
        w=barywei(theta);
        d2=2;
        dMDM_DDE = kron(D(2:end,:), eye(d2));
        PM=@(u,j,x)barint(theta,w,u(j:d2:end)',x);

    otherwise
        dMDM_DDE = kron(D(2:end,:), eye(1));
end

switch model
    case 'DL'  % Delay logistic equation
        dydt = [par(1)*u(1)*(1-u(end)/par(2)); dMDM_DDE*u];
        

    case 'MG'   % Mackey glass equation
        dydt = [par(1)*u(end)/(1+u(end)^par(3)) - par(2)*u(1); dMDM_DDE*u];

    case {'Rossler1','Rossler2'} % Rossler system with two cases: one delay case and two delays case
        dydt = [-u(2)-u(3)+par(1)*PM(u,1,-tau(1))+par(2)*u(end-2);...
            u(1)+par(3)*u(2);par(4)-par(5)*u(3)+u(1)*u(3);dMDM_DDE*u];

    case 'tau_3' % Two neuron system
        dydt = [-par(1)*u(1) + par(2)*tanh(PM(u,1,-tau(1))) + par(3)*tanh(PM(u,2,-tau(3)));...
            -par(1)*u(2) + par(2)*tanh(PM(u,2,-tau(1))) + par(4)*tanh(PM(u,1,-tau(2)));...
            dMDM_DDE*u];
        
    otherwise
        error('Invalid model type.');
end
end

