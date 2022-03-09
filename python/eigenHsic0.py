import numpy as np
from scipy import linalg as LA
import sortem

def eigenHsic0(X, Y, d):
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

	D, V = LA.eig(M)
	D = np.diag(D)
	D = D.real
	V = V.real

	V, D = sortem.sortem(V, D)

	D = D[:d,:d];
	V = V[:,:d];

	value = np.trace(D)/((m-1)**2)
	D = np.diag(D)
	return value, D, V