function f=hsca0(X,Y,d,lambda)



[n,m] = size(X);
[n1,m1] = size(Y);



if (m1 ~= m)
    error('Number of instances in X and Y must be equal');
end

A = ones(m,m);
H = eye(m,m) - (1/m)*A;



%
L = Y'*Y;
%L = (Y'*Y + 1).^2;
%L = exp(-0.1*distance(Y,Y));

M= X*H*L*H*X';
M1 = eye(size(M));
Lf = zeros(m,m);

S = [];
P = [];
for i=1:d
[V,D] = eig(M,M1);
D = real(D);
V = real(V);
[V,D] = sortem(V,D);
P(:,i) = V(:,1);
S(i) = D(1,1);

Lf = X'*P(:,1:i)*P(:,1:i)'*X;
M1 = (X*H*Lf*H*X' + lambda*eye(size(M1)));
end

f.value = sum(S);
f.D = diag(S);
f.V = P;

