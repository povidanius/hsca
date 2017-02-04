function f=hsca1(X,Y,d,lambda)

%X = X'; Y = Y'; d = 4;

[n,m] = size(X);
[n1,m1] = size(Y);



if (m1 ~= m)
    error('Number of instances in X and Y must be equal');
end

A = ones(m,m);
H = eye(m,m) - (1/m)*A;

L = Y'*Y;
%L = (Y'*Y + 1).^2;
%L = exp(-1*distance(Y,Y));

%M= X*H*L*H*X';
LL = L - diag(diag(L));
M = X*(LL + (A*LL*A - sum(sum(LL))*eye(m,m))/((m-1)*(m-2)) - (LL*A + A*LL - 2*diag(diag(LL*A)))/(m-2))*X'; 
Lf = zeros(m,m);
M1 = eye(size(M));
S = [];
P = [];
for i=1:d
    [V,D] = eig(M,M1);
D = real(D);
V = real(V);
[V,D] = sortem(V,D);
P(:,i) = V(:,1);
S(i) = D(1,1);
%Lf = Lf + X'*P(:,i)*P(:,i)'*X;
%M1 = (X*H*Lf*H*X' + 0.000001*eye(size(M1)));
%Lf = gram(X'*P(:,1:i),X'*P(:,1:i),'linear');
Lf = X'*P(:,1:i)*P(:,1:i)'*X;

LLf = Lf - diag(diag(Lf));

M1 = X*(LLf + (A*LLf*A - sum(sum(LLf))*eye(m,m))/((m-1)*(m-2)) - (LLf*A + A*LLf - 2*diag(diag(LLf*A)))/(m-2))*X' + lambda*eye(size(M1)); 
end

f.value = sum(S);
f.D = diag(S);
f.V = P;

