function Theta_all = poolData(X_data, polyorder,usesine)

    n = size(X_data,1);
    nVars=size(X_data,2);
%     disp(size(X_data));
%     disp(size(nVars));
    ind = 1;
    % poly order 0
    Theta_all(:,ind) = ones(n,1);
    ind = ind+1;
% disp(size(Theta_all));
    % poly order 1
    for i=1:nVars
        Theta_all(:,ind) = X_data(:,i);
        ind = ind+1;
    end

    if(polyorder>=2)
        % poly order 2
        for i=1:nVars
            for j=i:nVars
                Theta_all(:,ind) = X_data(:,i).*X_data(:,j);
                ind = ind+1;
            end
        end
    end

    if(polyorder>=3)
        % poly order 3
        for i=1:nVars
            for j=i:nVars
                for k=j:nVars
                    Theta_all(:,ind) = X_data(:,i).*X_data(:,j).*X_data(:,k);
                    ind = ind+1;
                end
            end
        end
    end

    

    if(polyorder>=4)
        % poly order 4
        for i=1:nVars
            for j=i:nVars
                for k=j:nVars
                    for l=k:nVars
                        Theta_all(:,ind) = X_data(:,i).*X_data(:,j).*X_data(:,k).*X_data(:,l);
                        ind = ind+1;
                    end
                end
            end
        end
    end

    if(polyorder>=5)
        % poly order 5
        for i=1:nVars
            for j=i:nVars
                for k=j:nVars
                    for l=k:nVars
                        for m=l:nVars
                            Theta_all(:,ind) = X_data(:,i).*X_data(:,j).*X_data(:,k).*X_data(:,l).*X_data(:,m);
                            ind = ind+1;
                        end
                    end
                end
            end
        end
    end
if usesine
    trig_terms = [];
    for k = 1:10
        trig_terms = [trig_terms, tanh(k * X_data)];
    end
    Theta_all = [Theta_all, trig_terms];
end

end
