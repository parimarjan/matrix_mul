#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

#define SEED 1234

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

/* used very often because struct matrix's data member is float ** while in the contiguous case we
 * use float *
 */
#define contig_data(A) ((float *) A->m)

#define NUM_PARALLEL_THREADS 4

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
        free(mat->m);
    } 
    free(mat);
}

/* does not allocate contiguous arrays. Using calloc so I dont have to set stuff to 0 later. */
struct matrix *alloc_matrix(int rows, int cols, enum alloc_type type) { 
     struct matrix *mat = malloc(sizeof(struct matrix));
     if (type == naive) {
         float **C = (float **)calloc(rows, sizeof(float*));
         int r;
         for (r = 0; r < rows; ++r) {
          C[r] = (float *)calloc(cols, sizeof(float));
         }
         mat->m = C;
     } else if (type == ac) {
        float **m;
        int i;
        m = (float**)calloc(rows, sizeof(float*));
        m[0] = (float*)calloc(rows * cols, sizeof(float));
        for (i = 1; i < rows; ++i)
            m[i] = m[i-1] + cols; 
        mat->m = m;
     } else if (type == contig) {
         mat->m= (float **) calloc(rows * cols, sizeof(float)) ;
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
     if (type == contig) {
        for (int i = 0; i < (rows * cols); i++) {
            /* annoying to recast... */
            contig_data(mat)[i] = (rand()+0.0)/(RAND_MAX+0.0);
        }
     } else {
         for (r = 0; r < rows; ++r) {
             for (c = 0; c < cols; ++c) {
                mat->m[r][c] = (rand()+0.0)/(RAND_MAX+0.0);
              }
         }
     }
     return mat;
}

/* TODO: switch order of j and k loops */
void mult_naive(struct matrix *C,
	     struct matrix *A,
         struct matrix *B)
{
     int i, j, k;
     int r;
     
     for (i = 0; i < A->rows; i++) {
	  for(j = 0; j < B->cols; j++) {
	       for(k = 0; k < A->cols; k++) {
		    C->m[i][j] += A->m[i][k] * B->m[k][j];
	       }
	    } 
     }
}

void mult_naive_kj(struct matrix *C,
	     struct matrix *A,
         struct matrix *B)
{
     int i, j, k;
     int r;
     
     for (i = 0; i < A->rows; i++) {
       for(k = 0; k < A->cols; k++) {
          for(j = 0; j < B->cols; j++) {
		    C->m[i][j] += A->m[i][k] * B->m[k][j];
	       }
	    } 
     }
}


/* attractivechaos stuff */
void ac_mat_mul0(struct matrix *C,
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
}

void ac_mat_mul1(struct matrix *C,
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
}

void ac_mat_mul1_kj(struct matrix *C,
                           struct matrix *A,
                           struct matrix *B)
{
    int i, j, k, n_b_rows = A->cols;
    struct matrix *bT;
    bT = mat_transpose(B->rows, B->cols, B);
    for (i = 0; i < A->rows; ++i) {
        const float *ai = A->m[i];
        float *mi = C->m[i];
        for (k = 0; k < A->cols; ++k)
            for (j = 0; j < B->cols; ++j) {
                mi[j] += ai[k] * bT->m[j][k];
        }
    }
    free_matrix(bT);
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
void ac_mat_mul2(struct matrix *C, struct matrix *A, struct matrix *B)
{
    int i, j;
    struct matrix *bT;
    bT = mat_transpose(B->rows, B->cols, B);
    for (i = 0; i < A->rows; ++i)
        for (j = 0; j < B->cols; ++j) 
            C->m[i][j] = sdot_sse(A->cols, A->m[i], bT->m[j]);
    free_matrix(bT);
}

void ac_mat_mul7(struct matrix *C, struct matrix *A, struct matrix *B)
{
    int x = 16;
	int i, j, ii, jj;
	struct matrix *bT = mat_transpose(B->rows, B->cols, B);
	for (i = 0; i < A->rows; i += x) {
		for (j = 0; j < B->cols; j += x) {
			int je = B->cols < j + x? B->cols : j + x;
			int ie = A->rows < i + x? A->rows : i + x;
			for (ii = i; ii < ie; ++ii)
				for (jj = j; jj < je; ++jj)
					C->m[ii][jj] += sdot_sse(A->cols, A->m[ii], bT->m[jj]);
		}
	}
	free_matrix(bT);
}
#endif

#ifdef HAVE_CBLAS
#include <cblas.h>

void ac_mat_mul6(struct matrix *C, struct matrix *A, struct matrix *A)
{
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A->rows, B->cols, A->cols, 1.0f, A->m[0], A->rows, B->m[0], B->rows, 0.0f, C->m[0], C->rows);
}
#endif


/*********************************************
 * end of the stuff from attractive chaos repo 
 *********************************************
*/

/* contiguous stuff */
void contig_blocked(struct matrix *C, struct matrix *A, struct matrix *B) {
    int block_size = A->rows / 8;
    /*int block_size = 16;*/
    for (int kk = 0; kk < C->rows; kk += block_size) {
        for (int ii = 0; ii < A->rows; ii += block_size) {
            for (int jj = 0; jj < B->cols; jj += block_size) {
                for (int k = kk; k < min(kk + block_size, C->rows); k++) {
                    for (int i = ii; i < min(ii + block_size, A->rows); i++) {
                        for (int j = jj; j < min(jj + block_size, B->cols); j++) {
                            contig_data(C)[i*C->rows + j] += contig_data(A)[i*A->rows + k] * contig_data(B)[k*B->cols + j];
                        }
                    }
                }
            }
        }
    }
}

/* TODO: Switch the order of j and k loops */
void contig_naive(struct matrix *C, struct matrix *A, struct matrix *B) {
    for (int i = 0; i < A->rows; i++) {
        for (int k = 0; k < C->rows; k++) {
            for (int j = 0; j < B->cols; j++) {
                contig_data(C)[i*C->rows + j] +=  contig_data(A)[i*A->rows + k] * contig_data(B)[k*A->rows + j];
            }
        }
    }
}

void get_thread_range(int total_work, int num_threads, int tid, int *start, int *end) {
    /*dividing up grace bucket indices based on tid. The only overlap between*/
    /*threads will be reading field values - but that should not be a prob. */
  *start = (total_work / num_threads) * tid;
  *end = *start + (total_work / num_threads);
  if (*end > total_work || tid == num_threads - 1) {
      *end = total_work;
  }
}

void thread_workload_contig_naive(struct matrix *C, struct matrix *A, struct matrix *B, int tid) {
    int start, end;
    get_thread_range(C->rows, NUM_PARALLEL_THREADS, tid, &start, &end);

    for (int i = start; i < end; i++) {
        for (int k = 0; k < C->rows; k++) {
            for (int j = 0; j < B->rows; j++) {
                contig_data(C)[i*C->rows + j] += contig_data(A)[i*A->rows + k] * contig_data(B)[k*B->rows + j];
            }
        }
    }
}

/* Parallelizing just the outer loop */
void contig_naive_par(struct matrix *C, struct matrix *A, struct matrix *B) {
#pragma omp parallel for
    for (int i = 0; i < NUM_PARALLEL_THREADS; i++) {
      thread_workload_contig_naive(C, A, B, i);
    }
}


