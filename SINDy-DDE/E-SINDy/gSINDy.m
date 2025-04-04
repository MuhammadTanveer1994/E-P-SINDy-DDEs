function dx = gSINDy(t, y, Z, ahat, polyorder,usesine, model,bestpar3)


    switch model
        case 'MG'
            X_data = [y(:)', Z(:)', 1./(1 + Z(:).^bestpar3)];
        otherwise
            X_data = [y(:)',Z(:)'];
    end

    Theta = poolData(X_data, polyorder,usesine);

    dx =(Theta * ahat)';
end
