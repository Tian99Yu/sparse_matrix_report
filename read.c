#include <stdio.h>
#include <stdlib.h>
#include "mmio.h"

int main(){
    MM_typecode matcode;
    FILE *f;
    f = fopen("./TSOPF_RS_b678_c2/TSOPF_RS_b678_c2.mtx", "r");
    mm_read_banner(f, &matcode);
    printf("haha %d", mm_is_sparse(matcode));
}