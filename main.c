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
        if (Lx[Lp[j]] != 1){
            printf("this diagonal value is wrong, val: %f, col %d, row %d", Lx[Lp[j]], j, Lp[j]);
        }
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
        //check if x[j] is 0 to save time
        if (x[j]==0){continue;}
        x[j] /= Lx[Lp[j]];
        for (p = Lp[j] + 1; p < Lp[j + 1]; p++)
        {
            x[Li[p]] -= Lx[p] * x[j];
        }
    }
    return (1);
}


int lsolve_improve_2(int n, int *Lp, int *Li, double *Lx, double *x)
{
    int p, j;
    double x_tmp;
    if (!Lp || !Li || !x)
        return (0);
    /* check inputs */
    for (j = 0; j < n; j++)
    {
        //check if x[j] is 0 to save time
        if (x[j]==0){continue;}
        x[j] /= Lx[Lp[j]];
        x_tmp = x[j];
        for (p = Lp[j] + 1; p < Lp[j + 1]; p++)
        {
          x[Li[p]] -= Lx[p] * x_tmp; }
    }
    return (1);
}



int lsolve_DFS_traversal(int n, int *Lp, int *Li, double *Lx, double *x)
{
    int p, j;
    if (!Lp || !Li || !x)
        return (0);
    /* check inputs */
    int * visited = malloc(n * sizeof(int));
    memset(visited, 0, sizeof(int) * n);
    int * stack = malloc(n * sizeof(int));
    int * non_zeros = malloc(n * sizeof(int));
    int stack_size, non_zeros_size;
    stack_size = 0;
    non_zeros_size=0;
    //init the stack to be non-zero elements of b
    for (int i=0; i<n; i++){
        if(x[i]!=0){
            stack[stack_size] = i;
            stack_size++;
        }
    }
    int cur_element, cur_col;
    //apply DFS on the matrix to find all adjacent nodes
    while(stack_size != 0){
        //pop the element
        cur_element = stack[stack_size -1];
        stack_size--;
        if (visited[cur_element]==1) continue;
        //push it into the X array
        non_zeros[non_zeros_size] = cur_element;
        non_zeros_size++;
        //mark it as visited
        visited[cur_element] = 1;
        //push all the node's non-visited children to the stack
        for (int i=Lp[cur_element]+1; i<Lp[cur_element +1]; i++){
            if (visited[Li[i]] == 0){                       
                stack[stack_size] = Li[i];
                stack_size++;
            }
        }
    }

    heapSort(non_zeros, non_zeros_size);
    for(int i=0; i<non_zeros_size; i++)
    {
        j=non_zeros[i];
        x[j] /= Lx[Lp[j]];
        for (p = Lp[j] + 1; p<Lp[j+1] ; p++)
        {
            x[Li[p]] -= Lx[p] * x[j];
        }
    }

    free(visited);
    free(stack);
    free(non_zeros);
    return (1);
}


int lsolve_improve_omp(int n, int *Lp, int *Li, double *Lx, double *x)
{
    int p, j;
    if (!Lp || !Li || !x)
        return (0);
    /* check inputs */
    for (j = 0; j < n; j++)
    {
        //check if x[j] is 0 to save time
        if (x[j]==0){continue;}
        x[j] /= Lx[Lp[j]];
        
        #pragma omp parallel private(p) num_threads(2)
        #pragma omp for
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

int verification(Matrix * mtx, double* b, double* answer){
    int nz, dim;
    nz = mtx->nz;
    dim = mtx->dim;
    double result[dim];
    memset(result, 0, sizeof(double) * dim);
    for (int i=0;i<dim; i++){
        if (answer[i] > 100000) printf("blow value, %f", answer[i]);
    }
    for(int col=0;col<dim; col++){
        for(int i=mtx->Lp[col]; i < mtx->Lp[col+1];i++){
            int row;
            double val;
            row = mtx->Li[i];
            val = mtx->Lx[i];

        if (answer[col] > 100000) printf("blow answer, %f", answer[col]);
        if (val > 100000) printf("blow value, %f", val);
        if (answer[col] * val > 100000) printf("blow multiplication %f",val * answer[col] );
            result[row] += answer[col] * val;
        }
    } 
    for (int i=0; i<dim;i++){
        if (abs(b[i] - result[i]) > 0.0001){
            printf("b: %f, result: %f, iteration %d\n", b[i], result[i], i);
            return 0;
        }
    }
    return 1;

}



int get_time(int (*solver_pt)(int, int*, int*, double*, double*), Matrix* m, double * solution, double* time, double * verification_b){
    //read matrices
    //get dimension
    struct timespec time_start, time_finish;
    clock_gettime(CLOCK_MONOTONIC, &time_start);
    (*solver_pt)(m->dim, m->Lp, m->Li, m->Lx, solution);
    clock_gettime(CLOCK_MONOTONIC, &time_finish);

    struct timespec time_diff = difftimespec(time_finish, time_start);
    * time = timespec_to_msec(time_diff);

    int validate_result;
    validate_result = verification(m, verification_b, solution);
    return validate_result;
}



int main(){

    // double *solution1, *solution2, *solution3,* verification_b;
    // Matrix * m1, *debug_m;
    // double t1, t2, t3;
    // m1 = read_matrix("matrices/debug_2/matrix.mtx");
    // read_b("matrices/debug_2/b.mtx", &solution1);
    // // lsolve_DFS_traversal(m1->dim, m1->Lp, m1->Li, m1->Lx, solution1);

    // read_b("matrices/debug_2/b.mtx", &verification_b);
    // int r3 = get_time(&lsolve_DFS_traversal, m1, solution1, &t3, verification_b);
    // for (int i=0; i<9; i++){
    //     printf("%f\n", solution1[i]);

    // }
    // printf("time %f, verification %d\n", t3, r3);


    //==================================================
    double *solution1, *solution2, *solution3,* verification_b;
    Matrix * m1, *debug_m;
    double t1, t2, t3;
    m1 = read_matrix("matrices/TSOPF_RS_b678_c2/TSOPF_RS_b678_c2.mtx");
    read_b("matrices/TSOPF_RS_b678_c2/b_for_TSOPF_RS_b678_c2_b.mtx", &solution1);
    read_b("matrices/TSOPF_RS_b678_c2/b_for_TSOPF_RS_b678_c2_b.mtx", &solution2);
    read_b("matrices/TSOPF_RS_b678_c2/b_for_TSOPF_RS_b678_c2_b.mtx", &solution3);
    read_b("matrices/TSOPF_RS_b678_c2/b_for_TSOPF_RS_b678_c2_b.mtx", &verification_b);

    int r1 = get_time(&lsolve, m1, solution1, &t1, verification_b);
    int r2 = get_time(&lsolve_improve_2, m1, solution2, &t2, verification_b);
    int r3 = get_time(&lsolve_DFS_traversal, m1, solution3, &t3, verification_b);
    printf("time %f, verification %d\n", t1, r1);
    printf("time %f, verification %d\n", t2, r2);
    printf("time %f, verification %d\n", t3, r3);
    exit(0);

}