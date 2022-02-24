import numpy as np
from scipy import linalg as LA
import sortem

def eigenHsic1(X, Y, d):
	n = X.shape[0]
	m = X.shape[1]
	n1 = Y.shape[0]
	m1 = Y.shape[1]


	if m != m1:
		print('Number of instances in X and Y must be equal')
		return

	A = np.ones((m, m))
	# H = np.eye(m) - A/m

	L = Y.getH() * Y
	LL = L - np.diag(np.diag(L))
	M = X * (LL + (A * LL * A - np.sum(np.sum(LL)) * np.eye(m))/((m-1)*(m-2)) - (LL * A + A * LL - 2 * np.diag(np.diag(LL * A)))/(m-2)) * X.getH()

	D, V = LA.eig(M)
	D = np.diag(D)
	D = D.real
	V = V.real
	V, D = sortem.sortem(V, D)

	D = D[:d,:d];
	V = V[:,:d];

	value = np.trace(D)/((m-3)*m)
	D = np.diag(D)
	return value, D, V