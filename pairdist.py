import numpy as np


pairdist = np.genfromtxt("PairwiseDistance.csv", delimiter=",")

print(pairdist[:6])

for i in range(1, pairdist.shape[0]):
    for j in range(1, pairdist.shape[1]):
        if i < j:
            pairdist[j, i] = pairdist[i, j]

print(pairdist[1:])

np.savetxt("pairdist.csv", pairdist[1:], delimiter=",")
