import numpy as np

def solver(L, b):
    x = b
    for j in range(L.shape[0]):
        i=j+1
        while i < L.shape[0]:
            if L[i][j] != 0:
                x[i] = x[i] - L[i][j] * x[j]
            i += 1
    return x


if __name__ == "__main__":
    L = np.array([[1, 0], [2, 3]])
    b = np.array([3,3])
    print(solver(L, b))