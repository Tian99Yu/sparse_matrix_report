#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "mmio.h"
#include "read.h"
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
/**
 * @brief improved version of lsolve
 * 
 * @param n 
 * @param Lp 
 * @param Li 
 * @param Lx 
 * @param x 
 * @return int 
 */
int lsolve_improve_1(int n, int *Lp, int *Li, double *Lx, double *x)
{
    int p, j;
    if (!Lp || !Li || !x)
        return (0);
    /* check inputs */
    for (j = 0; j < n; j++)
    {
        //check if x[j] is 0 to save time
        if (x[j] == 0)
        {
            continue;
        }
        x[j] /= Lx[Lp[j]];
        for (p = Lp[j] + 1; p < Lp[j + 1]; p++)
        {
            x[Li[p]] -= Lx[p] * x[j];
        }
    }
    return (1);
}

/**
 * @brief improved virsion 2 of lsolve, implemented the DFS traversal to find
 * the related dependencies
 * 
 * @param n 
 * @param Lp 
 * @param Li 
 * @param Lx 
 * @param x 
 * @return int 
 */
int lsolve_DFS_traversal(int n, int *Lp, int *Li, double *Lx, double *x)
{
    int p, j;
    if (!Lp || !Li || !x)
        return (0);
    /* check inputs */
    int *visited = malloc(n * sizeof(int));
    memset(visited, 0, sizeof(int) * n);
    int *stack = malloc(n * sizeof(int));
    int *non_zeros = malloc(n * sizeof(int));
    int stack_size, non_zeros_size;
    stack_size = 0;
    non_zeros_size = 0;
    //init the stack to be non-zero elements of b
    for (int i = 0; i < n; i++)
    {
        if (x[i] != 0)
        {
            stack[stack_size] = i;
            stack_size++;
        }
    }
    int cur_element;
    //apply DFS on the matrix to find all adjacent nodes
    while (stack_size != 0)
    {
        //pop the element
        cur_element = stack[stack_size - 1];
        stack_size--;
        if (visited[cur_element] == 1)
            continue;
        //push it into the X array
        non_zeros[non_zeros_size] = cur_element;
        non_zeros_size++;
        //mark it as visited
        visited[cur_element] = 1;
        //push all the node's non-visited children to the stack
        for (int i = Lp[cur_element] + 1; i < Lp[cur_element + 1]; i++)
        {
            if (visited[Li[i]] == 0)
            {
                stack[stack_size] = Li[i];
                stack_size++;
            }
        }
    }

    heapSort(non_zeros, non_zeros_size);
    for (int i = 0; i < non_zeros_size; i++)
    {
        j = non_zeros[i];
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
 * @brief Function used to verify the result from the lsolve function. 
 * It calculates a vector y=Lx and compare y with the original b
 * 
 * @param mtx The matrix L
 * @param b The vector b
 * @param answer The result from lsolve function
 * @return int 1 if the answer is correct and 0 otherwise
 */
int verification(Matrix *mtx, double *b, double *answer)
{
    int dim;
    dim = mtx->dim;
    //result is the multiplication result (y where y=Lx), I will
    //compare the vector result with  the original b later
    double result[dim];
    //process calculating the result vector using L*x
    memset(result, 0, sizeof(double) * dim);
    for (int col = 0; col < dim; col++)
    {
        for (int i = mtx->Lp[col]; i < mtx->Lp[col + 1]; i++)
        {
            int row;
            double val;
            row = mtx->Li[i];
            val = mtx->Lx[i];
            result[row] += answer[col] * val;
        }
    }
    //compare result with vector b
    int correct = 1;
    for (int i = 0; i < dim; i++)
    {
        if (abs(b[i] - result[i]) > 0.0001)
        {
            fprintf(stderr, "diff %f, b %f, r %f, i %d\n", b[i]- result[i], b[i], result[i], i);
            correct = 0;
        }
    }
    return correct;
}
/**
 * @brief Function used to calculate the time spent on the lsolve function
 * 
 * @param solver_pt The functin pointer, point to the lsolve function
 * @param m The matrix (L)
 * @param solution The vector, (b)
 * @param time The variable to record the time it takes
 * @param verification_b The same vector b, but used for verification
 * @return int The verification result, 1 for correct answer and 0 for wrong answer
 */
int get_time(int (*solver_pt)(int, int *, int *, double *, double *), Matrix *m, double *solution, double *time, double *verification_b)
{
    //read matrices
    //get dimension
    struct timespec time_start, time_finish;
    clock_gettime(CLOCK_MONOTONIC, &time_start);
    (*solver_pt)(m->dim, m->Lp, m->Li, m->Lx, solution);
    clock_gettime(CLOCK_MONOTONIC, &time_finish);

    struct timespec time_diff = difftimespec(time_finish, time_start);
    *time = timespec_to_msec(time_diff);

    int validate_result;
    validate_result = verification(m, verification_b, solution);
    return validate_result;
}

/**
 * @brief Function used to calculate the speedup
 * 
 * @param baseline 
 * @param cur_time 
 * @return double 
 */
double get_speedup(double baseline, double cur_time)
{
    return baseline / cur_time;
}


void print_m(Matrix *m)
{
    FILE *f = fopen("matrix.txt", "w");
    fprintf(f, "CCSMatrix\n");
    fprintf(f, "num_row: %d\n", m->dim);
    fprintf(f, "num_col: %d\n", m->dim);
    fprintf(f, "num_val: %d\n", m->nz);
    fprintf(f, "column pointer: \n");
    for (int i = 0; i < m->dim+1; i++)
    {
        fprintf(f, " %d ", m->Lp[i]);
    }
    fprintf(f, "\n");
    fprintf(f, "row index: \n");
    for (int i=0; i<m->nz; i++){
        fprintf(f, " %d ", m->Li[i]);
    }
}


int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        fprintf(stderr, "Usage: ./main.out **matrix dir** **b vector dir**");
        return 1;
    }
    char *mtx_dir = argv[1];
    char *b_dir = argv[2];

    double *solution1, *solution2, *solution3, *verification_b;
    double su1, su2, su3;
    Matrix *m1;

    double t1, t2, t3;
    m1 = read_matrix(mtx_dir);
    read_b(b_dir, &solution1);
    read_b(b_dir, &solution2);
    read_b(b_dir, &solution3);
    read_b(b_dir, &verification_b);

    int r1 = get_time(&lsolve, m1, solution1, &t1, verification_b);
    // int r2 = get_time(&lsolve_improve_1, m1, solution2, &t2, verification_b);
    // int r3 = get_time(&lsolve_DFS_traversal, m1, solution3, &t3, verification_b);
    su1 = get_speedup(t1, t1);
    // su2 = get_speedup(t1, t2);
    // su3 = get_speedup(t1, t3);
    printf("time %f, speed up %f, verification %d\n", t1, su1, r1);
    // printf("time %f, speed up %f, verification %d\n", t2, su2, r2);
    // printf("time %f, speed up %f, verification %d\n", t3, su3, r3);
    return 0;
}