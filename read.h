#ifndef READ_H
#define READ_H

typedef struct
{
    int * Li, *Lp;
    double * Lx;
    int dim, nz;
} Matrix;
Matrix* read_matrix(char *dir);
void read_b(char * dir, double** b);
int get_dim(char *dir);
void heapSort(int* arr, int size);

#endif