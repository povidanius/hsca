import numpy as np
import hcsa0
import hcsa1
import eigenHsic0
import eigenHsic1
import sortem
from scipy import linalg as LA

def load_data(filename):
	with open(filename, "r") as file:
	    return np.matrix([line.split() for line in file], dtype=float)


X = load_data("../data.txt")
Y = load_data("../y.txt")

lamda = 1e-5
numFeatures = 7

value, D, V = hcsa0.hcsa0(X.getH(), Y.getH(), numFeatures, lamda)
Xp = X * V
print(Xp)

value, D, V = hcsa1.hcsa1(X.getH(), Y.getH(), numFeatures, lamda)
Xp = X * V
print(Xp)

value, D, V = eigenHsic0.eigenHsic0(X.getH(), Y.getH(), numFeatures)
Xp = X * V
print(Xp)

value, D, V = eigenHsic1.eigenHsic1(X.getH(), Y.getH(), numFeatures)
Xp = X * V
print(Xp)
