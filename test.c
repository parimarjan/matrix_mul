#include"matrix.h"
#include"stdio.h"
#include<stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>

/* TODO: Support non-square matrices */
void print_matrix(float **m, int rows, int cols) {
     int r,c;	  
     printf("[ ");
     for (r = 0; r < rows; r++) {
	  printf("[ ");
	  for (c = 0; c < cols; c++) {
	       printf("%f", m[r][c]);
	       if (c != cols -1 ) {
		    printf(", ");
	       }
	  }
	  printf(" ]\n");
     }
     printf("]\n"); 
}

long compute_sum(struct matrix *mat) {
    float sum = 0;

    if (mat->type == ac || mat->type == naive) {
        for (int i = 0; i < mat->rows; ++i) {
            for (int j = 0; j < mat->cols; ++j) {
                sum += mat->m[i][j];
            }
        }
    } else {
        for (int i = 0; i < (mat->rows * mat->cols); i++) {
            /* annoying to recast... */
            sum += ((float *) mat->m)[i];
        }
    }
    return sum;
}

/* Runs an arbitrary multiply implementation and times it. 
 * @f: matrix multiplication function.
 * @n: size of matrix (assuming square for simplicity)
 * @block: block size
 */
void run_multiply(void (*f) (struct matrix *, struct matrix *, struct matrix*), int n, char *name, enum alloc_type type) 
{
     struct matrix *a = rand_matrix(n,n,type);
     struct matrix *b = rand_matrix(n,n,type);
     struct matrix *c = alloc_matrix(n,n,type);
    
     /* Time only the multiply function call */
     struct timeval start, end, diff;
     gettimeofday(&start, 0);
     f(c,a,b);

     gettimeofday(&end, 0);
     timersub(&end, &start, &diff);
     /* sanity check with computing sum of matrix */
     long sum = compute_sum(c);
     printf("%s took time: %ld.%06ld (result=%ld)\n",
        name, (long) diff.tv_sec, (long) diff.tv_usec, sum);

     /* Cleanup */
     free_matrix(a);
     free_matrix(b);
     free_matrix(c);
}

int main(int argc, char *argv[]) {
    
    /* Technically, need to support non-square matrices, but all the results we care about are only
     * square matrices, so just accept n as a parameter */
    int n = 128;
    int ch;

    while ((ch = getopt(argc, argv, "n:b:")) != -1) {
        switch (ch) {
            case 'n':
                n = atof(optarg);
                break;
            case '?':
            default:
                fprintf(stderr, "invalid options");
                exit(1);
        }
    }
    /* Need a smarter way to choose block size */ 
    printf("n = %d\n", n);
    /* run it with different optimization schemes */
    run_multiply(mult_naive, n, "C-Naive", naive);
    /*run_multiply(mult_naive_kj, n, "C-Naive-kj", naive);*/
    /*run_multiply(mult_parallel_kj, n, "C-Parallel-kj", naive);*/

    /* attractive chaos based implementation */
    /*run_multiply(ac_mat_mul0, n, "ac_naive", ac);*/
    /*run_multiply(ac_mat_mul1, n, "ac_transposed", ac);*/
    /*run_multiply(ac_mat_mul1_kj, n, "ac_transposed_kj", ac);*/

#ifdef __SSE__
    /*run_multiply(ac_mat_mul2, n, "ac_sse", ac);*/
    /*run_multiply(ac_mat_mul7, n, "ac_sse_blocked", ac);*/
#endif

/*#ifdef HAVE_CBLAS*/
    /*run_multiply(ac_mat_mul6, n, "ac_blas", ac);*/
/*#endif*/
    
    /* contiguous arrays */
    /*run_multiply(contig_naive_kj, n, "contig_naive", contig);*/
    /*run_multiply(contig_blocked, n, "contig_blocked", contig);*/
}
