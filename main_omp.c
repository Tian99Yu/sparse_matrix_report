#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "mmio.h"
#include "queue.h"
#include "read.h"
#include "time_util.h"

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
    int cur_element, cur_col;
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
        for (int i = Lp[cur_element]; i < Lp[cur_element + 1]; i++)
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

int lsolve_level_omp(int n, int *Lp, int *Li, double *Lx, double *x)
{
    //===============analysis phase ==================================

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
    for (int i = 0; i < n; i++)
    {
        if (x[i] != 0 && num_parent[i] != 1)
        {
            //there is a subtree top node
            restart_loop = 1;
            x[i] = 0;
        }
    }
    //Here is  the place to run DFS a second time for that purpose
    if (restart_loop)
    {
        destroyQueue(queue1);
        queue1 = createQueue(n);
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

    //==================compute phase===================================
    int index = 0;
    int arr[n];
    for (int i=0; i<level_pt_size; i++){
        //for each level
        int cur_upper = level_pt[i];
        int cur_size = cur_upper - index;
        for (int j=0;j<cur_size; j++){
            arr[j] = level[index + j];
        }
        //first sort the node in that level
        heapSort(arr, cur_size);
        //then process the children
        #pragma omp parallel default(shared) num_threads(16)
        #pragma omp for 
        for(int child=0; child<cur_size; child++)
        {
            int j=arr[child];
            x[j] /= Lx[Lp[j]];
            for (int p = Lp[j] + 1; p<Lp[j+1] ; p++)
            {
                x[Li[p]] -= Lx[p] * x[j];
            }
        }
        index=cur_upper;

    }
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
        if (x[j] == 0)
        {
            continue;
        }
        x[j] /= Lx[Lp[j]];

#pragma omp parallel default(shared) private(p) num_threads(16)
#pragma omp for nowait
        for (p = Lp[j] + 1; p < Lp[j + 1]; p++)
        {
            x[Li[p]] -= Lx[p] * x[j];
        }
    }
    return (1);
}

int verification(Matrix *mtx, double *b, double *answer)
{
    int nz, dim;
    nz = mtx->nz;
    dim = mtx->dim;
    double result[dim];
    memset(result, 0, sizeof(double) * dim);
    // for (int i=0;i<dim; i++){
    //     if (answer[i] > 100000) printf("blow value, %f", answer[i]);
    // }
    for (int col = 0; col < dim; col++)
    {
        for (int i = mtx->Lp[col]; i < mtx->Lp[col + 1]; i++)
        {
            int row;
            double val;
            row = mtx->Li[i];
            val = mtx->Lx[i];

            // if (answer[col] > 100000) printf("blow answer, %f", answer[col]);
            // if (val > 100000) printf("blow value, %f", val);
            // if (answer[col] * val > 100000) printf("blow multiplication %f",val * answer[col] );
            result[row] += answer[col] * val;
        }
    }
    FILE *f = fopen("./log.txt", "w");

    for (int i = 0; i < dim; i++)
    {
        if (abs(b[i] - result[i]) > 0.0001)
        {
            fprintf(f, "b: %f, result: %f, iteration %d\n", b[i], result[i], i);
            // return 0;
        }
    }
    return 1;
}

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

int main()
{

    double *solution1, *solution2, *solution3, *verification_b;
    Matrix *m1, *debug_m;
    double t1, t2, t3;
    // m1 = read_matrix("matrices/debug_2/matrix.mtx");
    // read_b("matrices/debug_2/b.mtx", &solution1);
    // read_b("matrices/debug_2/b.mtx", &verification_b);

    m1 = read_matrix("matrices/TSOPF_RS_b678_c2/TSOPF_RS_b678_c2.mtx");
    read_b("matrices/TSOPF_RS_b678_c2/b_for_TSOPF_RS_b678_c2_b.mtx", &solution1);
    read_b("matrices/TSOPF_RS_b678_c2/b_for_TSOPF_RS_b678_c2_b.mtx", &verification_b);
    int r3 = get_time(&lsolve_level_omp, m1, solution1, &t3, verification_b);
    printf("time %f, verification %d\n", t3, r3);

    //==================================================

    // double *solution1, *solution2, *solution3,* verification_b;
    // Matrix * m1, *debug_m;
    // double t1, t2, t3;
    // m1 = read_matrix("matrices/TSOPF_RS_b678_c2/TSOPF_RS_b678_c2.mtx");
    // read_b("matrices/TSOPF_RS_b678_c2/b_for_TSOPF_RS_b678_c2_b.mtx", &solution1);
    // read_b("matrices/TSOPF_RS_b678_c2/b_for_TSOPF_RS_b678_c2_b.mtx", &solution2);
    // read_b("matrices/TSOPF_RS_b678_c2/b_for_TSOPF_RS_b678_c2_b.mtx", &solution3);
    // read_b("matrices/TSOPF_RS_b678_c2/b_for_TSOPF_RS_b678_c2_b.mtx", &verification_b);

    // int r1 = get_time(&lsolve_improve_omp, m1, solution1, &t1, verification_b);
    // // int r2 = get_time(&lsolve_improve_1, m1, solution2, &t2, verification_b);
    // // int r3 = get_time(&lsolve_DFS_traversal, m1, solution3, &t3, verification_b);
    // printf("time %f, verification %d\n", t1, r1);
    // // printf("time %f, verification %d\n", t2, r2);
    // // printf("time %f, verification %d\n", t3, r3);
    // exit(0);
}