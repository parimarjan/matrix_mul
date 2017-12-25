#include"stdio.h"
#include<stdlib.h>

enum alloc_type {
    naive,
    contig,
    ac,
};

struct matrix {
    float **m;
    int rows;
    int cols;
    enum alloc_type type;
};

void free_matrix(struct matrix *);
struct matrix *alloc_matrix(int, int, enum alloc_type);
struct matrix *rand_matrix(int, int, enum alloc_type);

/* different matrix multiplication methods */
struct matrix *mult_naive(struct matrix *, struct matrix *, struct matrix *);


/* attractive chaos ones */
struct matrix *ac_mat_mul0(struct matrix *, struct matrix *, struct matrix *);
struct matrix *ac_mat_mul1(struct matrix *, struct matrix *, struct matrix *);
struct matrix *ac_mat_mul2(struct matrix *, struct matrix *, struct matrix *);
//struct matrix *ac_mat_mul3(struct matrix *, struct matrix *, struct matrix *);

struct matrix *ac_mat_mul6(struct matrix *, struct matrix *, struct matrix *);

