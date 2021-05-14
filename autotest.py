import subprocess
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.linalg
import math
def execute_command(command):
    # print(">>>>> Executing command {}".format(command), flush=True)
    process = subprocess.Popen(command, stdout = subprocess.PIPE,
                               stderr = subprocess.STDOUT, shell=True,
                               universal_newlines = True)
    return_code = process.wait()
    output = process.stdout.read()

    if return_code == 1:
        print("failed to execute command = ", command)
        print(output)
        # exit()

    return output

#generate matrix
def gen_matrix(dir, dim, sparse_rate):
    header = """%%MatrixMarket matrix coordinate real general
%-------------------------------------------------------------------------------
% UF Sparse Matrix Collection, Tim Davis
% http://www.cise.ufl.edu/research/sparse/matrices/Norris/torso1
% name: Norris/torso1
% [S.Norris, Univ Auckland. finite diff/boundary elem,  2D model of torso]
% id: 896
% date: 2003
% author: S. Norris
% ed: T. Davis
% fields: title A name id date author ed kind
% kind: 2D/3D problem
%-----------------------------------------------------------
"""
    rows = np.random.randint(0, dim, math.ceil(dim**2 * sparse_rate))
    cals = np.random.randint(0, dim, math.ceil(dim**2 * sparse_rate))
    vals = np.random.randint(1, 50, math.ceil(dim**2 * sparse_rate))
    mtx = np.zeros((dim,dim))
    with open(dir, "w") as f:
        for i in range(math.ceil(dim**2 * sparse_rate)):
            index = [rows[i], cals[i]]
            index.sort()
            c, r = index
            assert(r>=c)
            mtx[r,c] = vals[i]
        for i in range(dim):
            for j in range(dim):
                if i==j: mtx[i, j]=1
        nz=5
        for i in range(dim):
            for j in range(dim):
                if mtx[i,j] !=0:
                    nz+=1
        f.write(header)
        f.write("{} {} {}\n".format(dim, dim,nz ))
        for j in range(dim):
            for i in range(dim):
                if mtx[i,j] != 0:
                    f.write("{} {} {}\n".format(i+1,j+1, mtx[i,j]))
    return mtx


def gen_b(dir, dim, num_element):
    header = """%%MatrixMarket matrix coordinate real general
%-------------------------------------------------------------------------------
% UF Sparse Matrix Collection, Tim Davis
% http://www.cise.ufl.edu/research/sparse/matrices/Norris/torso1
% name: Norris/torso1
% [S.Norris, Univ Auckland. finite diff/boundary elem,  2D model of torso]
% id: 896
% date: 2003
% author: S. Norris
% ed: T. Davis
% fields: title A name id date author ed kind
% kind: 2D/3D problem
%-----------------------------------------------------------
"""
    index = np.random.randint(0, dim, num_element)
    index.sort()
    vals = np.random.randint(0, 100, num_element)
    b = np.zeros((dim,1))
    with open(dir, "w") as f:
        f.write(header)
        f.write("{} 1 {}\n".format(dim, num_element))
        for i in range(num_element):
            f.write("{} 1 {}\n".format(index[i]+1, vals[i]))
            b[index[i],0] = vals[i]
    return b

def parse_row(row):
    row = row.split(" ")
    return [int(row[0]), float(row[1])]
def validate_function(p, dim=20, sparse_rate=0.1):
    L = gen_matrix("./matrix.mtx", dim, sparse_rate)
    if p:
        print(L)
    b = gen_b("./b.mtx", dim, 3)
    if p:
        print(b)
    execute_command("rm ./*.txt")
    execute_command("make")
    ret = execute_command("./main.out ./matrix.mtx ./b.mtx")
    x = scipy.linalg.solve_triangular(L, b, lower=True, unit_diagonal=True)
    x = x.reshape(-1)
    x_cout = np.zeros(dim)
    with open("./torso_x_coutput.txt", "r") as f:
        l = f.readline()
        while l:
            index, val = parse_row(l)
            x_cout[index] = val
            l = f.readline()
    if p:
        print(ret)
        print(x)
        print(x_cout)
    for i in range(x.shape[0]):
        if x[i] != x_cout[i]:
            print("x: {} x_cout: {}".format(x[i], x_cout[i]))
    return np.sum(x_cout-x), np.sum(L@x_cout - b)


def validate_looping(num_loop):
    for i in range(num_loop):
        x_s, b_s = validate_function(False)
        print(x_s, b_s)
        if x_s != 0 or b_s != 0:
            return False
    return True


if __name__ == "__main__":
    np.set_printoptions(suppress=True)
    # print(validate_function(False, 1000))
    print(validate_looping(10000))