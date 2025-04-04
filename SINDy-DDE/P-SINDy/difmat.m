function [D, x] = difmat(a, b, M)
        % CHEB pseudospectral differentiation matrix on Chebyshev nodes.
        % [D,x]=CHEB(a,b,M) returns the pseudospectral differentiation matrix D
    
        if M == 0
            x = 1;
            D = 0;
            return
        end
        x = ((b - a) * cos(pi * (0:M)' / M) + b + a) / 2;
        % x = cos(pi*(0:M)/M)';
        c = [2; ones(M-1, 1); 2].*(-1).^(0:M)';
        X = repmat(x, 1, M+1);
        dX = X - X';
        D = (c * (1./c)')./(dX + (eye(M+1)));
        D = D - diag(sum(D'));
        
end


