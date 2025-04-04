function v=barywei(x)
%BARYWEI barycentric weights.
%  v=barywei(x) returns the weights v for barycentric interpolation at
%  the nodes x according to [1].
%
%   REFERENCES:
%   [1] J.P. Berrut and L.N. Trefethen,"Barycentric Lagrange
%       interpolation",SIAM Rev. 46(3):501-517,2004.

n=length(x); %number of nodes
dx=ones(n,1);
v=dx;
for m=2:n
    for i=1:m-1
        dx(i)=x(m)-x(i);
        v(i)=-dx(i)*v(i);
    end
    v(m)=prod(dx(1:m));
end
v=1./v;
v=v/max(v); %normalization to 1