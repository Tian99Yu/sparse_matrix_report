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
    //Li is the row index
    Li = (int *)malloc(nz * sizeof(int));
    //Lp is the column pointer
    Lp = (int *)malloc((N + 1) * sizeof(int));
    Lx = (double *)malloc(nz * sizeof(double));
    int cur_col;
    int row, col;
    double val;
    for (int i = 0; i < nz; i++)
    {

        fscanf(f, "%d %d %lg\n", &row, &col, &val);
        row--;
        col--;
        if (cur_col != col)
        {
            cur_col = col;
            (Lp)[cur_col] = i;
        }
        (Li)[i] = row;
        (Lx)[i] = val;
    }
    (Lp)[N] = nz;
    mtx->dim = M;
    mtx->Li = Li;
    mtx->Lp = Lp;
    mtx->Lx = Lx;
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
        row--;
        col--;
        (*b)[row] = val;
    }
}

// int main()
// {
//     // int *Li, *Lp;
//     // double *Lx;
//     // read_matrix("matrices/TSOPF_RS_b678_c2/TSOPF_RS_b678_c2.mtx", &Li, &Lp, &Lx);
//     // for (int i = 0; i < 5; i++)
//     // {
//     //     printf("line %d: row index %d, col pointer %d, val %f\n",i, Li[i], Lp[i], Lx[i]);
//     // }

//     // printf("last 2 position of col pointer %d, %d",Lp[35696-1], Lp[35696]);

//     // double *b;
//     // read_b("./matrices/TSOPF_RS_b678_c2/b_for_TSOPF_RS_b678_c2_b.mtx", &b);
//     // printf("random location %f, 169 pos %f, 1695 %f", b[4], b[168], b[1694]);
// }