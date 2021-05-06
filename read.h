#ifndef READ_H
#define READ_H

void read_matrix(char *dir, int **Li, int **Lp, double **Lx);
void read_b(char * dir, double** b);
int get_dim(char *dir);

#endif