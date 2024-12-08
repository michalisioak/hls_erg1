import numpy as np

A = np.array([
    [0.333333, 0.288675, 0.258199, 0, 0],
    [0.288675, 0.25, 0.223607, 0.25, 0],
    [0.258199, 0.223607, 0.2, 0.223607, 0.258199],
    [0, 0.25, 0.223607, 0.25, 0.288675],
    [0, 0, 0.258199, 0.288675, 0.333333],
])
H_in = np.loadtxt("feature_matrix.txt", delimiter=",", dtype=int)
weights = np.loadtxt("weights.txt", delimiter=",", dtype=int)
H_out = np.maximum(0,A @ H_in @ weights)

print("python_H_out:")
for i in range(H_out.shape[0]):
    for j in range(H_out.shape[1]):
        print(f"{H_out[i][j]:{2}.2e}", end=' ')
    print() 


