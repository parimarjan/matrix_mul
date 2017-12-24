#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

/* Frees the matrix */
void free_matrix(struct matrix *mat) {
    if (mat->type == naive) { 
        int r;
        for (r = 0; r < mat->rows; ++r) {
            free(mat->m[r]);
        }
        free(mat->m);
        free(mat);
    }
}

/* does not allocate contiguous arrays */
struct matrix *alloc_matrix_naive(int rows, int cols) {
     float **C = (float **)malloc(sizeof(float*)*rows);
     int r;
     for (r = 0; r < rows; ++r) {
      C[r] = (float *)malloc(sizeof(float)*cols);
     }
     struct matrix *mat = malloc(sizeof(struct matrix));
     mat->m = C;
     mat->rows = rows;
     mat->cols = cols;
     mat->type = naive;
     return mat;
}

struct matrix *rand_matrix(int rows, int cols) {
     int r;
     int c;
     struct matrix *mat = alloc_matrix_naive(rows, cols);
     for (r = 0; r < rows; ++r) {
	  for (c = 0; c < cols; ++c) {
	       mat->m[r][c] = (rand()+0.0)/(RAND_MAX+0.0);
	  }
     }
     return mat;
}

/* @C: output matrix.
 */
struct matrix *mult_naive(struct matrix *C,
	     struct matrix *A,
         struct matrix *B)
{
     int i, j, k;
     int r;
     
     for (i = 0; i < A->rows; i++) {
	  for(j = 0; j < B->cols; j++) {
	       C->m[i][j] = 0;
	       for(k = 0; k < A->cols; k++) {
		    C->m[i][j] += A->m[i][k] * B->m[k][j];
	       }
	    } 
     }
     return C;
}
