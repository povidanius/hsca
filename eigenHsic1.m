function f=eigenHsic1(X,Y,d)



[n,m] = size(X);
[n1,m1] = size(Y);

if (m1 ~= m)
    error('Number of instances in X and Y must be equal');
end

A = ones(m,m);

L = Y'*Y;
%L = (Y'*Y + 1).^2;
%L = exp(-1*distance(Y,Y));



LL = L - diag(diag(L));

M = X*(LL + (A*LL*A - sum(sum(LL))*eye(m,m))/((m-1)*(m-2)) - (LL*A + A*LL - 2*diag(diag(LL*A)))/(m-2))*X'; 
[V,D] = eig(M);
V=real(V);
D=real(D); 


[V,D] = sortem(V,D); 

%sum(diag(D) > 0)

D = D(1:d,1:d);
V = V(:,1:d);



f.value = trace(D)/(m*(m-3));
f.D = diag(D);
f.V = V;

