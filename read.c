#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mmio.h"
#include "read.h"
Matrix *read_matrix(char *dir)
{

    int *Li, *Lp;
    double *Lx;
    Matrix *mtx = malloc(sizeof(Matrix));

    MM_typecode matcode;
    FILE *f;
    f = fopen(dir, "r");
    // error checking
    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }

    int M, N, nz;
    mm_read_mtx_crd_size(f, &M, &N, &nz);
    //Li is the row index, there is an extra M since I might need to add 1 to the diagonal
    Li = (int *)malloc((nz + M) * sizeof(int));
    //Lp is the column pointer
    Lp = (int *)malloc((N + 1) * sizeof(int));
    memset(Lp, -1, sizeof(int) * (N+1));
    Lx = (double *)malloc((nz + M) * sizeof(double));
    int cur_col;
    cur_col = -1;
    int row, col, total_nz;
    double val;
    total_nz = 0;
    for(int i=0; i<nz;i++){
        fscanf(f, "%d %d %lg\n", &row, &col, &val);
        row--;
        col--;
        if (row <= col){continue;}
        //if change to a new column
        if (cur_col != col){
            // case 1, no diagonal gap
            if (cur_col + 1 == col){
                //directly add a diagonal value
                cur_col ++;
                Lp[cur_col] = total_nz;
                Li[total_nz] = col;
                Lx[total_nz] = 1;
                total_nz ++;
            }else{
                //case2 diagonal gap, then fill all col's diagonals
                while(cur_col < col){
                    cur_col++;
                    Lp[cur_col] = total_nz;
                    Li[total_nz] = cur_col;
                    Lx[total_nz] = 1;
                    total_nz ++;
                }
            }
        }
        Li[total_nz] = row;
        Lx[total_nz] = val;
        total_nz++;
    }
    //fill out the zero diagonals till the end
    while(cur_col<M-1){
        cur_col++;
        Lp[cur_col] = total_nz;
        Li[total_nz] = cur_col;
        Lx[total_nz] = 1;
        total_nz++;
    }
    Lp[cur_col+1] = total_nz+1;


    mtx->dim = M;
    mtx->Li = Li;
    mtx->Lp = Lp;
    mtx->Lx = Lx;
    mtx->nz = total_nz;
    return mtx;
}

int get_dim(char *dir)
{

    MM_typecode matcode;
    FILE *f;
    f = fopen(dir, "r");
    // error checking
    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }

    int M, N, nz;
    mm_read_mtx_crd_size(f, &M, &N, &nz);
    return M;
}

void read_b(char *dir, double **b)
{
    MM_typecode matcode;
    FILE *f;
    f = fopen(dir, "r");
    // error checking
    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }

    int M, N, nz;
    mm_read_mtx_crd_size(f, &M, &N, &nz);

    *b = malloc(M * sizeof(double));
    memset(*b, 0, sizeof(double) * M);

    int row, col;
    double val;
    for (int i = 0; i < nz; i++)
    {

        // fscanf(f, "%d %d %lg\n", &((*Li)[i]), &((*Lp)[i]), &((*Lx)[i]));
        fscanf(f, "%d %d %lg\n", &row, &col, &val);
        //only use the lower part of the matrix
        row--;
        col--;
        (*b)[row] = val;
    }
}

int temp;  

void heapify(int arr[], int size, int i)
{
    int largest = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;

    if (left < size && arr[left] > arr[largest])
        largest = left;

    if (right < size && arr[right] > arr[largest])
        largest = right;

    if (largest != i)
    {
        temp = arr[i];
        arr[i] = arr[largest];
        arr[largest] = temp;
        heapify(arr, size, largest);
    }
}

void heapSort(int arr[], int size)
{
    int i;
    for (i = size / 2 - 1; i >= 0; i--)
        heapify(arr, size, i);
    for (i = size - 1; i >= 0; i--)
    {
        temp = arr[0];
        arr[0] = arr[i];
        arr[i] = temp;
        heapify(arr, i, 0);
    }
}
// int main()
// {
//     // int *Li, *Lp;
//     // double *Lx;
//     Matrix * m;
//     m = read_matrix("matrices/trivial_eg/matrix_bug.mtx");
//     for (int i=0; i<3; i++){
//         printf("%f\n", m->Lx[i]);
//     }
//     printf("caonima");
//     int * Lp;
//     Lp = m->Lp;
//     // for (int i = 0; i < 5; i++)
//     // {
//     //     printf("line %d: row index %d, col pointer %d, val %f\n",i, Li[i], Lp[i], Lx[i]);
//     // }

//     // printf("last 2 position of col pointer %d, %d",Lp[35696-1], Lp[35696]);

//     // double *b;
//     // read_b("./matrices/TSOPF_RS_b678_c2/b_for_TSOPF_RS_b678_c2_b.mtx", &b);
//     // printf("random location %f, 169 pos %f, 1695 %f", b[4], b[168], b[1694]);
// }