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
    long sum = 0;
    for (int i = 0; i < mat->rows; ++i) {
        for (int j = 0; j < mat->cols; ++j) {
            if (mat->type == contig) {
                sum += ((float *)mat->m)[i*mat->rows + j]; 
            } else if (mat->type == naive || mat->type == ac) {
                sum += mat->m[i][j];
            }
        }
    }
    return sum;
}

/* Runs an arbitrary multiply implementation and times it. 
 * @f: matrix multiplication function.
 * @n: size of matrix (assuming square for simplicity)
 * @block: block size
 */
void run_multiply(struct matrix* (*f) (struct matrix *, struct matrix *, struct matrix*), int n, int block, char *name, enum alloc_type type) 
{
     struct matrix *a = rand_matrix(n,n, type);
     struct matrix *b = rand_matrix(n,n, type);
     struct matrix *c = alloc_matrix(n,n, type);
    
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
    int ch, b;

    while ((ch = getopt(argc, argv, "n:b:")) != -1) {
        switch (ch) {
            case 'n':
                n = atof(optarg);
                break;
            case 'b':
                b = atof(optarg);
                break;
            case '?':
            default:
                fprintf(stderr, "invalid options");
                exit(1);
        }
    }
    /* Need a smarter way to choose block size */
    b = n/2;
    
    printf("n = %d\n", n);
    /* run it with different optimization schemes */
    run_multiply(mult_naive, n, 0, "naive", naive);
    /* attractive chaos based implementation */
    run_multiply(ac_mat_mul0, n, 0, "ac0", ac);
    run_multiply(ac_mat_mul1, n, 0, "ac1", ac);
    run_multiply(ac_mat_mul2, n, 0, "ac2", ac);

#ifdef HAVE_CBLAS
    run_multiply(ac_mat_mul6, n, 0, "ac6", ac);
#endif

}
