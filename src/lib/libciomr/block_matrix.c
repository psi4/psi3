/*!
  \file block_matrix.c
*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

void bzero();

/*!
** block_matrix() : Allocates a contiguous block of memory for an array of
** doubles, allocates an array of pointers to the beginning of each row and
** returns the pointer to the first row pointer.  This allows transparent
** 2d-array style access, but keeps memory together such that the matrix 
** could be used in conjunction with FORTRAN matrix routines.
**
** Allocates memory for an n x m matrix and returns a pointer to the
** first row. 
**
** T. Daniel Crawford
** Sometime in 1994
**
** Based on init_matrix() from libciomr
**
*/
double ** block_matrix(int n,int m)
   {
    double **A=NULL;
    double *B=NULL;
    int i;

    if(!m || !n) return((double **) NULL);

    if ((A = (double **) malloc(n * sizeof(double *)))==NULL) {
         fprintf(stderr,"block_matrix: trouble allocating memory \n");
         fprintf(stderr,"n = %d\n",n);
         exit(1);
         }

    if ((B = (double *) malloc(m*n * sizeof(double)))==NULL) {
         fprintf(stderr,"block_matrix: trouble allocating memory \n");
         fprintf(stderr,"m = %d\n",m);
         exit(1);
         }

    bzero(B, m*n*sizeof(double));

    for (i = 0; i < n; i++) {
         A[i] = &(B[i*m]);
         }

    return(A);
   }

void free_block(double **array)
   {
     if(array == NULL) return;
      free(array[0]);
      free(array);
   }
