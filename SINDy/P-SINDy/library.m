function V=library(n,g)
%n=number of variables (including state and nonlinear terms)
%e.g., x1,...,xd,sin(x1),...,sin(xd),1/x1,...,1/xd,...
%g=polynomial degree
%each row of V contains the indices of the variables appearing as factors
%in the considered monomial, e.g., V(i,:)=[1 1 2 3 3 4] represents the monomial
%X1*X1*X2*X3*X3*X4=X1^2*X2*X3^2*X4
v=[1:n]'; %local library degree 1
V=[0;v];  %gloabl library, degree 0 and 1
len(1)=size(v,1); %length local library
for k=2:g %cycle over degree k
    ind0=1;
    for i=1:size(v,1) %for any element of the local library
        for j=v(i,end):n %add all elements from the last index to n
            w(ind0,:)=[v(i,:),j]; %new local library of degree k
            ind0=ind0+1; %increase local length
        end
    end
    v=w; %update local library to new degree k
    V(sum(len)+1,k)=0; %augment gloabl library with one zero column
    V=[V;v]; %update global library with new degree k
    len(k)=size(v,1); %number of added elements of degree k
    w(:,size(w,2)+1)=0; %augment next local library with one zero column
end
size(V,1); %check of number of monomials of degree up to g in n variables
nchoosek(n+g,g)