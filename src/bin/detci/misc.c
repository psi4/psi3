/*! \file 
    \ingroup (DETCI)
    \brief Enter brief description of file here 
*/
#include <stdio.h>

/*
** TRANSP_SIGMA(): Function adds the transpose (times a phase factor) of
**    a matrix to itself.
**
*/
void transp_sigma(double **a, int rows, int cols, int phase)
{

   int i,j;

   if (rows != cols) {
      printf("(transp_sigma): Error, rows != cols\n");
      printf("\trows = %d, cols = %d\n", rows, cols);
      return;
      }

   /* do lower triangle */
   if (phase == 1) {
      for (i=0; i<rows; i++) {
         for (j=0; j<=i; j++) {
            a[i][j] += a[j][i];
            }
         }
      }
   else if (phase == -1) {
      for (i=0; i<rows; i++) {
         for (j=0; j<=i; j++) {
            a[i][j] -= a[j][i];
            }
         }
      }

   /* fix upper triangle (could remove me later if don't use upper tri) */
   if (phase == 1) {
      for (i=0; i<rows; i++) {
         for (j=i; j<cols; j++) {
            a[i][j] = a[j][i];
            }
         }
      }
   else {
      for (i=0; i<rows; i++) {
         for (j=i; j<cols; j++) {
            a[i][j] = -a[j][i];
            }
         }
      }

}


/*
** SET_ROW_PTRS()
**
** This function sets the row pointers for a 2D matrix allocated as one
** contiguous block of memory
**
*/
void set_row_ptrs(int rows, int cols, double **matrix)
{
   int i;
   double *ptr;

   ptr = matrix[0];

   for (i=1; i<rows; i++) {
      matrix[i] = matrix[0] + i * cols;
      }

}


