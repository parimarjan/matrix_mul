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
void mult_naive(struct matrix *, struct matrix *, struct matrix *);
void mult_naive_kj(struct matrix *, struct matrix *, struct matrix *);
void mult_parallel_kj(struct matrix *, struct matrix *, struct matrix *);

/* attractive chaos ones */
void ac_mat_mul0(struct matrix *, struct matrix *, struct matrix *);

void ac_mat_mul1(struct matrix *, struct matrix *, struct matrix *);
void ac_mat_mul1_kj(struct matrix *, struct matrix *, struct matrix *);

void ac_mat_mul2(struct matrix *, struct matrix *, struct matrix *);

void ac_mat_mul6(struct matrix *, struct matrix *, struct matrix *);
void ac_mat_mul7(struct matrix *, struct matrix *, struct matrix *);

/* contiguous allocation from weld */
void contig_naive_kj(struct matrix *, struct matrix *, struct matrix *);
void contig_blocked(struct matrix *, struct matrix *, struct matrix *);
