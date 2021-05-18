#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "mmio.h"
#include "queue.h"
#include "read.h"
#include "time_util.h"

// #ifndef DEBUG
// #define DEBUG
// #endif


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
 * @brief A refactored version of the lsolve_level_omp. However, the result is not correct. I
 * am still debugging it. It is not shown in the report. 
 * @param n 
 * @param Lp 
 * @param Li 
 * @param Lx 
 * @param x 
 */
int lsolve_level_improved_omp(int n, int *Lp, int *Li, double *Lx, double *x)
{
    #ifdef TIMER
    struct timespec t_s, t_f, t_d;
    #endif
    struct Queue *queue = createQueue(n);
    //the array of level information for each related node
    int level[n];
    int max_level = -1;
    memset(level, 0, sizeof(int) * n);
    for (int i = 0; i < n; i++)
    {
        if (x[i] != 0)
        {
            enqueue(queue, i);
            level[i] = 1;
        }
    }
    #ifdef TIMER
    clock_gettime(CLOCK_MONOTONIC, &t_s);
    #endif
    while (!isEmpty(queue))
    {
        int v = dequeue(queue);
        for (int i = Lp[v] + 1; i < Lp[v + 1]; i++)
        {
            level[Li[i]] = level[Li[i]] > level[v] + 1 ? level[Li[i]] : level[v] + 1;
            max_level = max_level > level[Li[i]] ? max_level : level[Li[i]];
            enqueue(queue, Li[i]);
        }
    }
    #ifdef TIMER
    clock_gettime(CLOCK_MONOTONIC, &t_f);
    t_d = difftimespec(t_f, t_s);
    fprintf(stderr, "BFS time %f", timespec_to_msec(t_d));
    #endif
    int level_info[max_level][n];
    int level_pt[max_level];
    #ifdef DEBUG
    int shit = 0;
    for (int i = 0; i < n; i++)
    {
        if (level[i] > 0)
            shit++;
    }
    fprintf(stderr, "total node involved %d", shit);
    #endif
    memset(level_pt, 0, sizeof(int) * n);
    int cur_level;
    for (int i = 0; i < n; i++)
    {
        cur_level = level[i];
        level_info[cur_level][level_pt[cur_level]] = i;
        level_pt[cur_level]++;
    }
    #ifdef TIMER
    clock_gettime(CLOCK_MONOTONIC, &t_s);
    #endif
    for (int i = 1; i < max_level + 1; i++)
    {
        int *arr = level_info[i];
        int cur_size = level_pt[i];
        for (int child = 0; child < cur_size; child++)
        {
            int j = arr[child];
            x[j] /= Lx[Lp[j]];
            for (int p = Lp[j] + 1; p < Lp[j + 1]; p++)
            {
                x[Li[p]] -= Lx[p] * x[j];
            }
        }
    }
    #ifdef TIMER
    clock_gettime(CLOCK_MONOTONIC, &t_f);
    t_d = difftimespec(t_f, t_s);
    fprintf(stderr, "compute time, %f\n", timespec_to_msec(t_d));
    #endif
#ifdef DEBUG
    int i = 1;
    FILE *f = fopen("level_log.txt", "w");
    int xx = 0;
    while (xx < max_level)
    {
        fprintf(f, "\n level %d =======\n", i);
        heapSort(level_info[i], level_pt[i]);
        for (int j = 0; j < level_pt[i]; j++)
        {
            fprintf(f, " %d ", level_info[i][j]);
        }
        i++;
        xx++;
    }
    // exit(1);
    printf("max level: %d", max_level);
    for (int i = 0; i < n; i++)
    {
        fprintf(stderr, "%d: %d\n", i, level[i]);
    }
#endif
    return 1;
}
/**
 * @brief The improved function from the naive omp approach.
 * It first generate the level sets based on the topological order
 * And then parallelize the computation of each level
 * 
 * param meanings are the same as lsolve
 * @param n 
 * @param Lp 
 * @param Li 
 * @param Lx 
 * @param x 
 * @return int 
 */
int lsolve_level_omp(int n, int *Lp, int *Li, double *Lx, double *x)
{
    //===============analysis phase ==================================
#ifdef TIMER
    struct timespec time_start1, time_finish1;
    struct timespec time_start2, time_finish2;
    clock_gettime(CLOCK_MONOTONIC, &time_start1);
#endif
    if (!Lp || !Li || !x)
        return (0);
    //init variables
    //level is used to store the level information of the dependency graph
    int *level = malloc(n * sizeof(int));
    int *level_pt = malloc(n * sizeof(int));
    //variables used in modified DFS & BFS
    int *stack = malloc(n * sizeof(int));
    int *num_parent = malloc(n * sizeof(int));
    int stack_size = 0;
    int level_size = 0;
    int level_pt_size = 0;
    int cur_element;
    struct Queue *queue1 = createQueue(n);
    struct Queue *queue2 = createQueue(n);
    memset(num_parent, 0, sizeof(int) * n);
    //init level 0 with all non-zero elements in b
    for (int i = 0; i < n; i++)
    {
        if (x[i] != 0)
        {
            // level[level_size] = i;
            stack[stack_size] = i;
            enqueue(queue1, i);
            stack_size++;
            // level_size ++;
        }
    }
    level_pt[level_pt_size] = level_size;
    //do a customized DFS (might visit one node multiple times)
    //the purpose is to get the number of parents of each node that should be processed
    while (stack_size != 0)
    {
        //pop the element
        cur_element = stack[stack_size - 1];
        stack_size--;
        if (num_parent[cur_element] > 0)
        {
            num_parent[cur_element] += 1;
            continue;
        }
        for (int i = Lp[cur_element] + 1; i < Lp[cur_element + 1]; i++)
        {
            if (num_parent[Li[i]] == 0)
            {
                //debug
                // if (cur_element == Li[i]){fprintf(stderr,"same element pushed and poped");}
                //end debug
                stack[stack_size] = Li[i];
                stack_size++;
            }
            else
            {
                num_parent[Li[i]]++;
            }
        }
        num_parent[cur_element]++;
    }

    //However, there is a case when some non-zero element of the input "b" is actually a
    //children node of other non-zero element of input "b", thus I need to run DFS a second
    //time after getting rid of them
    int restart_loop = 0;
    int not_enqueue[n];
    memset(not_enqueue, 0, n*sizeof(int));
    for (int i=0; i<n; i++){
        if (x[i]!=0 && num_parent[i]!=1){
            not_enqueue[i] = 1;
            restart_loop =1;
        }
    }
    // for (int i = 0; i < n; i++)
    // {
    //     if (x[i] != 0 && num_parent[i] != 1)
    //     {
    //         //there is a subtree top node
    //         restart_loop = 1;
    //         x[i] = 0;
    //     }
    // }
    //Here is  the place to run DFS a second time for that purpose
    if (restart_loop)
    {
        destroyQueue(queue1);
        queue1 = createQueue(n);
        memset(num_parent, 0, sizeof(int) * n);
        //init level 0 with all non-zero elements in b
        for (int i = 0; i < n; i++)
        {
            if (x[i] != 0 && not_enqueue[i]==0)
            {
                // level[level_size] = i;
                stack[stack_size] = i;
                enqueue(queue1, i);
                stack_size++;
                // level_size ++;
            }
        }
        level_pt[level_pt_size] = level_size;
        //do a customized DFS (might visit one node multiple times)
        //the purpose is to get the number of parents of each node that should be processed
        while (stack_size != 0)
        {
            //pop the element
            cur_element = stack[stack_size - 1];
            stack_size--;
            if (num_parent[cur_element] > 0)
            {
                num_parent[cur_element] += 1;
                continue;
            }
            for (int i = Lp[cur_element] + 1; i < Lp[cur_element + 1]; i++)
            {
                if (num_parent[Li[i]] == 0)
                {
                    //debug
                    if (cur_element == Li[i])
                    {
                        fprintf(stderr, "same element pushed and poped");
                    }
                    //end debug
                    stack[stack_size] = Li[i];
                    stack_size++;
                }
                else
                {
                    num_parent[Li[i]]++;
                }
            }
            num_parent[cur_element]++;
        }
    }

    //do a customized BFS to generate the level array of the dependency graph directly using the
    //given matrix
    //Later, each level could be processed in parallel (after sorting the nodes in each level)
    int iteration_count = 0;
    while (!isEmpty(queue1) || !isEmpty(queue2))
    {
        if (isEmpty(queue2))
        {
            while (!isEmpty(queue1))
            {
                cur_element = dequeue(queue1);
                if (num_parent[cur_element] != 1)
                {
                    fprintf(stderr, "cur_element %d, iteration count %d", cur_element, iteration_count);
                }
                //all elements in the queue should have one parent dependency left only
                for (int i = Lp[cur_element] + 1; i < Lp[cur_element + 1]; i++)
                {
                    //there are two cases, only node with 1 parent dependency will be pushed to the other queue
                    if (num_parent[Li[i]] == 1)
                    {
                        enqueue(queue2, Li[i]);
                    }
                    else
                    {
                        num_parent[Li[i]]--;
                    }
                }
                num_parent[cur_element]--;
                //now the node's num_parent should be 0
                //push this element to the current level
                if (num_parent[cur_element] != 0)
                {
                    printf("non zero!!!!!!");
                    exit(1);
                }
                level[level_size] = cur_element;
                level_size++;
            }
            //all elements are dequeued, level finish
            level_pt[level_pt_size] = level_size;
            level_pt_size++;
            iteration_count++;
        }
        else
        {
            while (!isEmpty(queue2))
            {
                cur_element = dequeue(queue2);
                //all element in queue2 should have only 1 parent dependency
                for (int i = Lp[cur_element] + 1; i < Lp[cur_element + 1]; i++)
                {
                    //there are two cases, only node with 1 parent dependency will be pushed to the other queue
                    if (num_parent[Li[i]] == 1)
                    {
                        enqueue(queue1, Li[i]);
                    }
                    else
                    {
                        //otherwise, we substract on dependency
                        num_parent[Li[i]]--;
                    }
                }
                num_parent[cur_element]--;
                //now push the node to the current level
                level[level_size] = cur_element;
                level_size++;
            }
            level_pt[level_pt_size] = level_size;
            level_pt_size++;
        }
    }
#ifdef TIMER
    clock_gettime(CLOCK_MONOTONIC, &time_finish1);

    struct timespec time_diff1 = difftimespec(time_finish1, time_start1);
    double t1 = timespec_to_msec(time_diff1);
    fprintf(stderr, "analysis time, %f", t1);

    clock_gettime(CLOCK_MONOTONIC, &time_start2);
#endif
    //==================compute phase===================================
    int index = 0;
    int arr[n];
#ifdef DEBUG
    FILE *f = fopen("level_correct.txt", "w");
    for (int i = 0; i < level_pt_size; i++)
    {
        int cur_upper = level_pt[i];
        int cur_size = cur_upper - index;
        for (int j = 0; j < cur_size; j++)
        {
            arr[j] = level[index + j];
        }
        //first sort the node in that level
        heapSort(arr, cur_size);
        fprintf(f, "\n level %d ======\n", i + 1);
        for (int m = 0; m < cur_size; m++)
        {
            fprintf(f, " %d ", arr[m]);
        }
        index = cur_upper;
    }
    index=0;
    // exit(1);
#endif

    for (int i = 0; i < level_pt_size; i++)
    {
        //for each level
        int cur_upper = level_pt[i];
        int cur_size = cur_upper - index;
        for (int j = 0; j < cur_size; j++)
        {
            arr[j] = level[index + j];
        }
        // //first sort the node in that level
        // heapSort(arr, cur_size);
//then process the children
#pragma omp parallel default(shared)
#pragma omp for
        for (int child = 0; child < cur_size; child++)
        {
            int j = arr[child];
            //no need to add omp critical here as elements at the same level are independent 
            x[j] /= Lx[Lp[j]];
            for (int p = Lp[j] + 1; p < Lp[j + 1]; p++)
            {
                #pragma omp critical
                x[Li[p]] -= Lx[p] * x[j];
            }
        }
        index = cur_upper;
    }
#ifdef TIMER
    clock_gettime(CLOCK_MONOTONIC, &time_finish2);
    struct timespec time_diff2 = difftimespec(time_finish2, time_start2);
    double t2 = timespec_to_msec(time_diff2);
    fprintf(stderr, "compute time, %f", t2);
#endif
    return 1;
}
/**
 * @brief The naive parallelization of the lsolve function
 * 
 * @param n 
 * @param Lp 
 * @param Li 
 * @param Lx 
 * @param x 
 * @return int 
 */
int lsolve_improve_omp(int n, int *Lp, int *Li, double *Lx, double *x)
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

#pragma omp parallel default(shared) private(p)
#pragma omp for
        for (p = Lp[j] + 1; p < Lp[j + 1]; p++)
        {
            x[Li[p]] -= Lx[p] * x[j];
        }
    }
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
    #ifdef DEBUG
    FILE *f = fopen("torso_x_coutput.txt", "w");
    for (int i = 0; i < mtx->dim; i++)
    {
        if (answer[i] == 0)
            continue;
        fprintf(f, "%d %lf\n", i, answer[i]);
    }
    #endif
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
            #ifdef DEBUG
            fprintf(stderr, "diff %f, b %f, r %f, i %d\n", b[i]- result[i], b[i], result[i], i);
            #endif
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
int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        fprintf(stderr, "Usage: ./main.out **matrix dir** **b vector dir**");
        return 1;
    }
    char *mtx_dir = argv[1];
    char *b_dir = argv[2];

    double serial_baseline;

    double *solution1, *solution2, *solution3, *verification_b;
    double su1, su2;
    Matrix *m1;
    double t1, t2;
    m1 = read_matrix(mtx_dir);
    read_b(b_dir, &solution1);
    read_b(b_dir, &solution2);
    read_b(b_dir, &solution3);
    read_b(b_dir, &verification_b);
    int r1 = get_time(&lsolve_improve_omp, m1, solution1, &t1, verification_b);
    int r2 = get_time(&lsolve_level_omp, m1, solution2, &t2, verification_b);
    get_time(&lsolve, m1, solution3, &serial_baseline, verification_b); 
    su1 = get_speedup(serial_baseline, t1);
    su2 = get_speedup(serial_baseline, t2);
    printf("time %f, speed up % f, verification %d\n", t1, su1, r1);
    printf("time %f, speed up %f, verification %d\n", t2, su2, r2);
    return 0;
}