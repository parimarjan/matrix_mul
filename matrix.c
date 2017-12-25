#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

#define SEED 1234

/* From the attractive chaos repository.
 */
struct matrix *mat_transpose(int n_rows, int n_cols, struct matrix *a)
{
    int i, j;
    if (a->type == ac) {
        struct matrix *mat = alloc_matrix(n_cols, n_rows, ac);
        for (i = 0; i < n_rows; ++i)
            for (j = 0; j < n_cols; ++j)
                mat->m[j][i] = a->m[i][j];
        return mat;
    }
    return NULL;
}

/* Frees the matrix */
void free_matrix(struct matrix *mat) {
    if (mat->type == naive) { 
        int r;
        for (r = 0; r < mat->rows; ++r) {
            free(mat->m[r]);
        }
        free(mat->m);
    } else if (mat->type == ac) {
        free(mat->m[0]); 
        free(mat->m);
    } else if (mat->type == contig) {

    } 
    free(mat);
}

/* does not allocate contiguous arrays */
struct matrix *alloc_matrix(int rows, int cols, enum alloc_type type) {
     
     struct matrix *mat = malloc(sizeof(struct matrix));
     if (type == naive) {
         float **C = (float **)malloc(sizeof(float*)*rows);
         int r;
         for (r = 0; r < rows; ++r) {
          C[r] = (float *)malloc(sizeof(float)*cols);
         }
         mat->m = C;
     } else if (type == ac) {
        float **m;
        int i;
        m = (float**)malloc(rows * sizeof(float*));
        m[0] = (float*)calloc(rows * cols, sizeof(float));
        for (i = 1; i < rows; ++i)
            m[i] = m[i-1] + cols; 
        mat->m = m;
     } else if (type == contig) {
        /* TODO */
     }

     mat->rows = rows;
     mat->cols = cols;
     mat->type = type;
     return mat;
}

struct matrix *rand_matrix(int rows, int cols, enum alloc_type type) {
     int r;
     int c;
     struct matrix *mat = alloc_matrix(rows, cols, type);
     srand(SEED);
     for (r = 0; r < rows; ++r) {
         for (c = 0; c < cols; ++c) {
             if (type == naive || type == ac) {
                mat->m[r][c] = (rand()+0.0)/(RAND_MAX+0.0);
             }
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

/* attractivechaos implementation */
struct matrix *ac_mat_mul0(struct matrix *C,
                           struct matrix *A,
                           struct matrix *B)
{
	int i, j, k;
	for (i = 0; i < A->rows; ++i) {
		for (j = 0; j < B->cols; ++j) {
			float t = 0.0;
			for (k = 0; k < A->cols; ++k)
				t += A->m[i][k] * B->m[k][j];
			C->m[i][j] = t;
		}
	}
	return C;
}

struct matrix *ac_mat_mul1(struct matrix *C,
                           struct matrix *A,
                           struct matrix *B)
{
    int i, j, k, n_b_rows = A->cols;
    struct matrix *bT;
    bT = mat_transpose(B->rows, B->cols, B);
    for (i = 0; i < A->rows; ++i) {
        const float *ai = A->m[i];
        float *mi = C->m[i];
        for (j = 0; j < B->cols; ++j) {
            float t = 0.0f, *bTj = bT->m[j];
            for (k = 0; k < A->cols; ++k)
                t += ai[k] * bTj[k];
            mi[j] = t;
        }
    }
    free_matrix(bT);
    return C;
}

#ifdef __SSE__
#include <xmmintrin.h>

float sdot_sse(int n, const float *x, const float *y)
{
	int i, n8 = n>>3<<3;
	__m128 vs1, vs2;
	float s, t[4];
	vs1 = _mm_setzero_ps();
	vs2 = _mm_setzero_ps();
	for (i = 0; i < n8; i += 8) {
		__m128 vx1, vx2, vy1, vy2;
		vx1 = _mm_loadu_ps(&x[i]);
		vx2 = _mm_loadu_ps(&x[i+4]);
		vy1 = _mm_loadu_ps(&y[i]);
		vy2 = _mm_loadu_ps(&y[i+4]);
		vs1 = _mm_add_ps(vs1, _mm_mul_ps(vx1, vy1));
		vs2 = _mm_add_ps(vs2, _mm_mul_ps(vx2, vy2));
	}
	for (s = 0.0f; i < n; ++i) s += x[i] * y[i];
	_mm_storeu_ps(t, vs1);
	s += t[0] + t[1] + t[2] + t[3];
	_mm_storeu_ps(t, vs2);
	s += t[0] + t[1] + t[2] + t[3];
	return s;
}
#endif

#ifdef __SSE__
struct matrix *ac_mat_mul2(struct matrix *C, struct matrix *A, struct matrix *B)
{
    int i, j;
    struct matrix *bT;
    bT = mat_transpose(B->rows, B->cols, B);
    for (i = 0; i < A->rows; ++i)
        for (j = 0; j < B->cols; ++j)
            C->m[i][j] = sdot_sse(A->cols, A->m[i], bT->m[j]);
    free_matrix(bT);
    return C;
}
#endif

#ifdef HAVE_CBLAS
#include <cblas.h>

/*float **mat_mul5(int n_a_rows, int n_a_cols, float *const *a, int n_b_cols, float *const *b)*/
/*{*/
	/*int i, j, n_b_rows = n_a_cols;*/
	/*float **m, **bT;*/
	/*m = mat_init(n_a_rows, n_b_cols);*/
	/*bT = mat_transpose(n_b_rows, n_b_cols, b);*/
	/*for (i = 0; i < n_a_rows; ++i)*/
		/*for (j = 0; j < n_b_cols; ++j)*/
			/*m[i][j] = cblas_sdot(n_a_cols, a[i], 1, bT[j], 1);*/
	/*mat_destroy(bT);*/
	/*return m;*/
/*}*/

struct matrix *ac_mat_mul6(struct matrix *C, struct matrix *A, struct matrix *A)
{
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A->rows, B->cols, A->cols, 1.0f, A->m[0], A->rows, B->m[0], B->rows, 0.0f, C->m[0], C->rows);
	return m;
}
#endif


/*********************************************
 * end of the stuff from attractive chaos repo 
 *********************************************
*/

