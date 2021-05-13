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
        
if __name__ == "__main__":
    t_f = input("T or F?")
    if t_f == "T":
        mdir = "matrices/TSOPF_RS_b678_c2/TSOPF_RS_b678_c2.mtx"
        bdir = "matrices/TSOPF_RS_b678_c2/TSOPF_RS_b678_c2.mtx"
    elif t_f == "F":
        mdir = "matrices/torso1/torso1.mtx"
        bdir = "matrices/torso1/b_for_torso1.mtx"
    else:
        mdir="matrices/debug_3/matrix.mtx"
        bdir="matrices/debug_3/b.mtx"
    L = read_mtx(mdir)
    b = read_b(bdir)
    x = spsolve(L, b)
    x = x.reshape((-1,1))
    with open("./torso1_x.txt", "w") as f:
        for i in range(x.shape[0]):
            if x[i, 0] != 0:
                f.write("{} {}\n".format(i, x[i,0]))

    # print(np.sum(L@x -b))
    # scipy.io.mmwrite("./torso1_x.mtx", x.reshape((-1,1)))
