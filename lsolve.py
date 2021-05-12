import numpy as np
from scipy.sparse.linalg import spsolve
import scipy.io
from scipy.sparse import csc_matrix

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
        rows, cols, vals = [], [], []
        cur_col = -1
        for i in range(nz):
            line = f.readline()
            i, j, v = parse_row(line)
            if j != cur_col:
                while cur_col+1<j:
                    cur_col+=1
                    rows.append(cur_col)
                    cols.append(cur_col)
                    vals.append(1)
                cur_col+=1
                if i != j:
                    rows.append(cur_col)
                    cols.append(cur_col)
                    vals.append(1)
                
            rows.append(i)
            cols.append(j)
            vals.append(v)
        while cur_col < nr-1:
            cur_col+=1
            rows.append(cur_col)
            cols.append(cur_col)
            vals.append(1)
        return csc_matrix((np.array(vals), (np.array(rows), np.array(cols))), shape=(nr,nr))


def read_b(dir):
    with open(dir, "r") as f:
        line = f.readline()
        while line[0] == "%":
            line = f.readline()
        header = line.split(" ")
        header = [int(i) for i in header]
        nr, nc, nz = header[0], header[1], header[2]
        rows, cols, vals = [], [], []
        for i in range(nz):
            line = f.readline()
            i, j, v = parse_row(line)
            rows.append(i)
            cols.append(0)
            vals.append(v)
        return csc_matrix((np.array(vals), (np.array(rows), np.array(cols))), shape=(nr,1))
        


L = read_mtx("matrices/torso1/torso1.mtx")
b = read_b("matrices/torso1/b_for_torso1.mtx")
x = spsolve(L, b)
x = x.reshape((-1,1))
print(np.sum(L@x -b))
# print(x.shape, type(x))
# scipy.io.mmwrite("./torso1_x.mtx", x.reshape((-1,1)))
