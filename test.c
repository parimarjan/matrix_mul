#include"matrix.h"
#include"stdio.h"
#include<stdlib.h>
#include <string.h>
/*#include <stdio.h>*/
#include <unistd.h>

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
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

/* Runs an arbitrary multiply implementation and times it. 
 * @f: matrix multiplication function.
 * @n: size of matrix (assuming square for simplicity)
 * @block: block size
 */
void run_multiply(float ** (*f) (float **, float **, int, int, float **, int, int), 
                  int n, int block) 
{
     printf("in run multiply!\n");
     float **a = rand_matrix(n,n);
     float **b = rand_matrix(n,n);
     float **c = alloc_matrix(n,n);
     f(c,a,n,n,b,n,n);
     printf("f done!\n");
     print_matrix(c,n,n);
     /* Cleanup */
     /*free(a);*/
     /*free(b);*/
     /*free(c);*/
}

int main(int argc, char *argv[]) {
    
    /* Technically, need to support non-square matrices, but all the results we care about are only
     * square matrices, so just accept n as a parameter */

    int n = 128;
    /* Need a smarter way to choose blocl size */
    int b = n/2;

    int ch;
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
    
    printf("n = %d\n", n);

    run_multiply(mult_naive, n, b);
}
