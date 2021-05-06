#include <stdio.h>
#include <stdlib.h>
#include "mmio.h"
#include "read.h"
int main()
{
    MM_typecode matcode;
    FILE *f;
    f = fopen("matrices/TSOPF_RS_b678_c2/TSOPF_RS_b678_c2.mtx", "r");
    mm_read_banner(f, &matcode);

    int M, N, nz;
    mm_read_mtx_crd_size(f, &M, &N, &nz);
    printf("M %d, N %d, nz %d", M, N, nz);

    // scan the first row
    int row, col;
    double val;
    fscanf(f, "%d %d %lg\n", &row, &col, &val);
    printf("row %d, col %d, val %f", row, col, val);
}