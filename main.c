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
        for (int i=Lp[cur_element]; i<Lp[cur_element +1]; i++){
            if (visited[Li[i]] == 0){                       
                stack[stack_size] = Li[i];
                stack_size++;
            }
        }
    }



    // for (j = 0; j < n; j++)
    for(int i=0; i<non_zeros_size; i++)
    {
        j=non_zeros[i];
        x[j] /= Lx[Lp[j]];
        for (p = Lp[j] + 1; p < Lp[j + 1]; p++)
        {
            x[Li[p]] -= Lx[p] * x[j];
        }
    }

    free(visited);
    free(stack);
    free(non_zeros);
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


int get_time(int (*solver_pt)(int, int*, int*, double*, double*), Matrix* m,char * b_dir, double * solution, double* time){
    //read matrices
    double *b;
    read_b(b_dir, &b);
    //get dimension
    struct timespec time_start, time_finish;
    clock_gettime(CLOCK_MONOTONIC, &time_start);
    (*solver_pt)(m->dim, m->Lp, m->Li, m->Lx, b);
    clock_gettime(CLOCK_MONOTONIC, &time_finish);

    struct timespec time_diff = difftimespec(time_finish, time_start);
    * time = timespec_to_msec(time_diff);

    int validate_result;
    validate_result = validation(solution, b, m->dim);
    free(b);
    return validate_result;
}



int main(){
    double *solution;
    Matrix * m1;
    m1 = read_matrix("matrices/TSOPF_RS_b678_c2/TSOPF_RS_b678_c2.mtx");
    read_b("./matrices/TSOPF_RS_b678_c2/b_for_TSOPF_RS_b678_c2_b.mtx", &solution);
    //get dimension
    int dim = get_dim("matrices/TSOPF_RS_b678_c2/TSOPF_RS_b678_c2.mtx");
    lsolve(m1->dim, m1->Lp, m1->Li, m1->Lx, solution);
    double t1, t2, t3;
    int return_1, return_2, return_3;
    return_1 = get_time(&lsolve, m1, "matrices/TSOPF_RS_b678_c2/b_for_TSOPF_RS_b678_c2_b.mtx", solution, &t1);
    printf("the time takes for 1 is %f, validate result %d\n", t1, return_1);
    return_2 = get_time(&lsolve_improve_1, m1, "matrices/TSOPF_RS_b678_c2/b_for_TSOPF_RS_b678_c2_b.mtx", solution, &t2);
    printf("the time takes for 1 is %f, validate result %d\n", t2, return_2);

    return_3 = get_time(&lsolve_DFS_traversal, m1, "matrices/TSOPF_RS_b678_c2/b_for_TSOPF_RS_b678_c2_b.mtx", solution, &t3);
    printf("the time takes for 1 is %f, validate result %d\n", t3, return_3);
}