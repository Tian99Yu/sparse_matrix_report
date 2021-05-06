#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "read.h"
#include "mmio.h"
#include "time_util.h"
/*
* Lower triangular solver Lx=b
* L is stored in the compressed column storage format
* Inputs are:
* n : the matrix dimension
* Lp : the column pointer of L
* Li : the row index of L
* Lx : the values of L
* In/Out:
* x : the right hand-side b at start and the solution x at the end.
*/
int lsolve(int n, int *Lp, int *Li, double *Lx, double *x)
{
    int p, j;
    if (!Lp || !Li || !x)
        return (0);
    /* check inputs */
    for (j = 0; j < n; j++)
    {
        x[j] /= Lx[Lp[j]];
        for (p = Lp[j] + 1; p < Lp[j + 1]; p++)
        {
            x[Li[p]] -= Lx[p] * x[j];
        }
    }
    return (1);
}

int main(){

    //read matrices
    int * Li, *Lp;
    double * Lx, *b;
    read_matrix("matrices/TSOPF_RS_b678_c2/TSOPF_RS_b678_c2.mtx", &Li, &Lp, &Lx);
    read_b("./matrices/TSOPF_RS_b678_c2/b_for_TSOPF_RS_b678_c2_b.mtx", &b);
    //get dimension
    int dim = get_dim("matrices/TSOPF_RS_b678_c2/TSOPF_RS_b678_c2.mtx");

    struct timespec time_start, time_finish;
    clock_gettime(CLOCK_MONOTONIC, &time_start);
    lsolve(dim, Lp, Li, Lx, b);
    clock_gettime(CLOCK_MONOTONIC, &time_finish);
    struct timespec time_diff = difftimespec(time_finish, time_start);
    double time_msec = timespec_to_msec(time_diff);
    printf("The time is %f\n", time_msec);
}