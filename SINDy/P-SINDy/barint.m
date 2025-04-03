function ff=barint(x,b,f,xx)
%BARINT Barycentric interpolation
%  Compute the values ff of a function f on xx using the barycentric
%  interpolation formula with x interpolation nodes, b barycentric
%  weights and f values of the function on x according to [1].
%
%  ff=BARINT(x,b,f,xx)
%
%   REFERENCES:
%   [1] J.P. Berrut and L.N. Trefethen,"Barycentric Lagrange
%       interpolation",SIAM Rev. 46(3):501-517,2004.

n=length(x);

numer=zeros(size(xx));
denom=zeros(size(xx));
exact=zeros(size(xx));
for j=1:n
    tdiff=xx-x(j);
    temp=b(j)./tdiff;
    numer=numer+temp*f(j);
    denom=denom+temp;
    exact(tdiff==0)=j;
end
jj=find(exact);
ff=numer./denom;
ff(jj)=f(exact(jj));