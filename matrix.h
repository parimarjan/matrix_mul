#include"stdio.h"
#include<stdlib.h>

enum alloc_type {
    naive,
    contig,
    two_malloc
};

struct matrix {
    float **m;
    int rows;
    int cols;
    enum alloc_type type;
};

void free_matrix(struct matrix *);
struct matrix *alloc_matrix_naive(int, int);
struct matrix *rand_matrix(int, int);

/* different matrix multiplication methods */
struct matrix *mult_naive(struct matrix *, struct matrix *, struct matrix *);
