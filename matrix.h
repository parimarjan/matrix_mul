#include"stdio.h"
#include<stdlib.h>

struct matrix {
    float **m;
    int rows;
    int columns;
};

float **alloc_matrix(int, int);
float **rand_matrix(int, int);

/* different matrix multiplication methods */
float **mult_naive(float **, float **, int, int, float **, int, int);
