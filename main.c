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

int lsolve_improve_1(int n, int *Lp, int *Li, double *Lx, double *x)
{
    int p, j;
    if (!Lp || !Li || !x)
        return (0);
    /* check inputs */
    for (j = 0; j < n; j++)
    {
        x[j] /= Lx[Lp[j]];
        //check if x[j] is 0 to save time
        if (x[j]==0){continue;}
        for (p = Lp[j] + 1; p < Lp[j + 1]; p++)
        {
            x[Li[p]] -= Lx[p] * x[j];
        }
    }
    return (1);
}


/**
 * @brief validate the output of our improved functions
 * 
 * @param solution the baseline solution,  an array of int
 * @param answer the answer given by improved functions
 * @param size the size of the integer array
 * @return return 0 on success and 1 on failure, ie: the solution is not correct
 */
int validation(double *solution, double * answer, int size){
    for (int i=0; i<size; i++){
        if (solution[i] != answer[i]){
            return 1;
        }
    }
    return 0;
}


int get_time(int (*solver_pt)(int, int*, int*, double*, double*), char * mtx_dir,char * b_dir, double * solution, double* time){
    //read matrices
    int * Li, *Lp;
    double * Lx, *b;
    read_matrix(mtx_dir, &Li, &Lp, &Lx);
    read_b(b_dir, &b);
    //get dimension
    int dim = get_dim(mtx_dir);
    struct timespec time_start, time_finish;
    clock_gettime(CLOCK_MONOTONIC, &time_start);
    (*solver_pt)(dim, Lp, Li, Lx, b);
    clock_gettime(CLOCK_MONOTONIC, &time_finish);

    struct timespec time_diff = difftimespec(time_finish, time_start);
    * time = timespec_to_msec(time_diff);

    int validate_result;
    validate_result = validation(solution, b, dim);
    free(Li);
    free(Lp);
    free(Lx);
    free(b);
    return validate_result;
}



int main(){
    int * Li, *Lp;
    double * Lx, *solution;
    read_matrix("matrices/TSOPF_RS_b678_c2/TSOPF_RS_b678_c2.mtx", &Li, &Lp, &Lx);
    read_b("./matrices/TSOPF_RS_b678_c2/b_for_TSOPF_RS_b678_c2_b.mtx", &solution);
    //get dimension
    int dim = get_dim("matrices/TSOPF_RS_b678_c2/TSOPF_RS_b678_c2.mtx");
    lsolve(dim, Lp, Li, Lx, solution);
    free(Li);
    free(Lp);
    free(Lx);
    double t1;
    int return_1;
    return_1 = get_time(&lsolve, "matrices/TSOPF_RS_b678_c2/TSOPF_RS_b678_c2.mtx", "matrices/TSOPF_RS_b678_c2/b_for_TSOPF_RS_b678_c2_b.mtx", solution, &t1);
    printf("the time takes for 1 is %f", t1);
}