static char *rcsid = "$Id$";

#include "includes.h"

extern void ludcmp(double **, int, int *, double *);
extern void lubksb(double **, int, int *, double *);
extern double *init_array(int);

/* solves linear equations A * x = b */
/* on input a contains coefficients, and b contains */
/* known vectors.  in is the dimension of a(in*in) */
/* im is the number of b vectors */
/* det returns determinant of matrix a */

void flin(a,b,in,im,det)
   double **a, *b, *det;
   int in, im;

{
int i,j,k,*indx;

    indx = (int *) init_array(in);

   ludcmp(a,in,indx,det);

   for (i=0; i < in ; i++) *det *= a[i][i];

   for (j=0; j<im; j++)
      lubksb(a,in,indx,b+j*in);

   free(indx);
}

