import numpy as np
from scipy import linalg as LA
import sortem

def hcsa0(X, Y, d, lamda):
	n = X.shape[0]
	m = X.shape[1]
	n1 = Y.shape[0]
	m1 = Y.shape[1]


	if m != m1:
		print('Number of instances in X and Y must be equal')
		return

	A = np.ones((m, m))
	H = np.eye(m) - A/m

	L = Y.getH() * Y
	M = X * H * L * H * X.getH()
	M1 = np.eye(M.shape[0])
	Lf = np.zeros((m, m))
	P = np.matrix(np.zeros((n, d)))
	S = np.zeros(d)
	for i in range(d):
		D, V = LA.eig(M, M1)
		D = np.diag(D)
		D = D.real
		V = V.real
		V, D = sortem.sortem(V, D)
		P[:,i] = np.reshape(V[:,0], (n, 1))
		S[i] = D[0,0]
		Lf = X.getH() * P[:,0:i+1] * P[:,0:i+1].getH() * X
		M1 = X * H * Lf * H * X.getH() + lamda * np.eye(n)

	value = np.sum(S)
	D = np.diag(S)
	V = P
	return value, D, V