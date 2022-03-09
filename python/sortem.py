import numpy as np

def sortem(P, D):
	ind = np.argsort(-np.diag(D))
	P2 = P[:,ind]
	D2 = np.diag(np.diag(D)[ind])
	return P2, D2