import subprocess
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def execute_command(command):
    print(">>>>> Executing command {}".format(command), flush=True)
    process = subprocess.Popen(command, stdout = subprocess.PIPE,
                               stderr = subprocess.STDOUT, shell=True,
                               universal_newlines = True)
    return_code = process.wait()
    output = process.stdout.read()

    if return_code == 1:
        print("failed to execute command = ", command)
        print(output)
        exit()

    return output

dirs = [("matrices/torso1/torso1.mtx","matrices/torso1/b_for_torso1.mtx" ), ("matrices/TSOPF_RS_b678_c2/TSOPF_RS_b678_c2.mtx", "matrices/TSOPF_RS_b678_c2/b_for_TSOPF_RS_b678_c2_b.mtx")]

def parse_row(r):
    r = r.split(",")
    speedup = r[1].split(" ")
    speedup = [i for i in speedup if len(i)>0]
    return float(speedup[2])

def plot_serial():
    execute_command("make")
    acc = []
    for i in range(len(dirs)):
        ret = execute_command("./main.out {} {}".format(dirs[i][0], dirs[i][1]))
        rows = ret.split("\n")
        rows = [r.strip() for r in rows if len(r)>0]
        rows = [parse_row(r) for r in rows]
        acc.append(rows)
    index = ["base line", "tiny improvement", "DFS version"]
    df = pd.DataFrame({
        "Torso": acc[0],
        "TSOPF": acc[1]
    }, index=index)

    ax = df.plot.bar()
    plt.title("Serial Speedup Comparison")
    plt.ylabel("Speedup")
    plt.xticks(rotation = 15)
    for p in ax.patches:
        ax.annotate(str(p.get_height())[:4], (p.get_x()+0.1, p.get_height() + 0.005), fontsize=8)
    plt.savefig("plots/{}.pdf".format("serial_speedup"))

def plot_omp():
    execute_command("make omp")
    acc = []
    for i in range(len(dirs)):
        ret = execute_command("./main_omp.out {} {}".format(dirs[i][0], dirs[i][1]))
        rows = ret.split("\n")
        rows = [r.strip() for r in rows if len(r)>0]
        rows = [parse_row(r) for r in rows]
        acc.append(rows)
    index = ["tiny improvement", "level version"]
    df = pd.DataFrame({
        "Torso": acc[0],
        "TSOPF": acc[1]
    }, index=index)

    ax = df.plot.bar()
    plt.title("Parallel Speedup Comparison, #T=16")
    plt.ylabel("Speedup")
    plt.xticks(rotation = 15)
    for p in ax.patches:
        ax.annotate(str(p.get_height())[:4], (p.get_x()+0.1, p.get_height() + 0.005), fontsize=8)
    plt.savefig("plots/{}.pdf".format("parallel_speedup"))
def plot_omp_diff_thread():
    execute_command("make omp")
    naive, level = [], []
    index = [2,4,8,16]
    for i in index:
        ret = execute_command("OMP_NUM_THREADS={} ./main_omp.out {} {}".format(i,dirs[1][0], dirs[1][1]))
        rows = ret.split("\n")
        rows = [r.strip() for r in rows if len(r)>0]
        rows = [parse_row(r) for r in rows]
        naive.append(rows[0])
        level.append(rows[1])
    df = pd.DataFrame({
        "naive":naive,
        "level":level
    }, index=index)
    df.plot.line()
    plt.title("The influence of #T on speed up")
    plt.ylabel("speed up")
    plt.xlabel("#T")
    plt.savefig("plots/t_influence.pdf")




if __name__== "__main__":
    plot_serial()
    plot_omp()
    plot_omp_diff_thread()
