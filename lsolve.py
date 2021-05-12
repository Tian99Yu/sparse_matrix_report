import numpy as np
from scipy.linalg import solve_triangular
import scipy

def parse_row(r):
    r = r.split(" ")
    return [int(r[0])-1, int(r[1])-1, float(r[2])] 


def read_mtx(dir):
    with open(dir, "r") as f:
        line = f.readline()
        while line[0] == "%":
            line=f.readline()
        header = line.split(" ")
        header = [int(i) for i in header]
        nr, nc, nz = header[0], header[1], header[2]
        #init the matrix
        mtx = np.zeros((nr, nc), dtype="float32")
        for i in range(nz):
            line = f.readline()
            i, j, v = parse_row(line)
            if(i >= j):
                mtx[i, j] = v
        for i in range(nr):
            if mtx[i,i] == 0:
                mtx[i,i] = 1
        return mtx

def read_b(dir):
    with open(dir, "r") as f:
        line = f.readline()
        while line[0] == "%":
            line = f.readline()
        header = line.split(" ")
        header = [int(i) for i in header]
        nr, nc, nz = header[0], header[1], header[2]
        b = np.zeros(nr, dtype="float32")
        for i in range(nz):
            line = f.readline()
            i, j, v = parse_row(line)
            b[i] = v
        return b


L = read_mtx("matrices/torso1/torso1.mtx")
b = read_b("matrices/torso1/b_for_torso1.mtx")
result = solve_triangular(L, b, lower=True)
with open("result.txt", "w") as f:
    for i in len(result):
        if result[i] != 0:
            f.write("{} {}\n".format(i, result[i]))
        
