#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

/* Frees the matrix */
void clean_matrix(struct matrix *m) {

}

float** alloc_matrix(int rows, int cols) {
     float **C = (float **)malloc(sizeof(float*)*rows);
     int r;
     for (r = 0; r < rows; ++r) {
	  C[r] = (float *)malloc(sizeof(float)*cols);
     }
     struct matrix *m = malloc(sizeof(struct matrix));
     return C;
}

float** rand_matrix(int rows, int cols) {
     int r;
     int c;
     float **C = alloc_matrix(rows, cols);
     for (r = 0; r < rows; ++r) {
	  for (c = 0; c < cols; ++c) {
	       C[r][c] = (rand()+0.0)/(RAND_MAX+0.0);
	  }
     }
     return C;
}

/* @C: output matrix.
 */
float** mult_naive(float **C,
	     float **A,
	     int a_rows,
	     int a_cols,
	     float **B,
	     int b_rows,
	     int b_cols) {
     int i, j, k;
     int r;
     
     for (i = 0; i < a_rows; i++) {
	  for(j = 0; j < b_cols; j++) {
	       C[i][j] = 0;
	       for(k = 0; k < a_cols; k++) {
		    C[i][j] += A[i][k] * B[k][j];
	       }
	    } 
     }

     return C;
}
