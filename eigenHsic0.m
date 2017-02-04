function f=eigenHsic0(X,Y,d)

%X = X'; Y = Y'; d = 4;

[n,m] = size(X);
[n1,m1] = size(Y);




if (m1 ~= m)
    error('Number of instances in X and Y must be equa');
end

A = ones(m,m);
H = eye(m,m) - (1/m)*A;

L = Y'*Y;
%L = (Y'*Y + 1).^2;
%L = exp(-1*distance(Y,Y));

M= X*H*L*H*X';



[V,D] = eig(M);
D = real(D);
V = real(V);
[V,D] = sortem(V,D);
%diag(D)

D = D(1:d,1:d);
V = V(:,1:d);

f.value = trace(D)/((m-1)^2);
f.D = diag(D);
f.V = V;

