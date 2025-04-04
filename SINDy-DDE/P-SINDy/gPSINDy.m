function dx = gPSINDy(t,u,ahat,polyorder,usesine,model,D,bestpar3)
switch model
    case 'MG'
         par(3) = bestpar3;  
end

  
    switch model
        case 'MG'
            X_data = [u', 1./ (1 + u.^par(3))'];
        otherwise
            X_data = u';
    end
    Theta = poolData(X_data, polyorder,usesine);
    uu = (Theta* ahat)';
    switch model
        case {'Rossler1','Rossler2'}
            dMDM_DDE =kron(D(2:end,:), eye(3));
            dx = [uu(1:3,:); dMDM_DDE*u];
            case 'tau_3'
            dMDM_DDE =kron(D(2:end,:), eye(2));
            dx = [uu(1:2,:); dMDM_DDE*u];
    otherwise
            dMDM_DDE = kron(D(2:end,:), eye(1));
            dx = [uu(1,:); dMDM_DDE*u];
    end
     
end
