
/* $Log$
 * Revision 1.1  2000/02/04 22:50:48  evaleev
 * Initial revision
 *
/* Revision 1.1  1991/06/15 22:45:28  seidl
/* Initial revision
/* */

static char *rcsid = "$Id$";

#include "includes.h"

linear_eq(a,na,nb,det)
   double **a,*det;
   int na,nb;

   {
      int i,j,*indx;
      double *b;

      b = (double *) init_array(na);
      indx = (int *) init_array(na);

      ludcmp(a,na,indx,det);

      for (i=0; i < na ; i++) *det *= a[i][i];

      for (i=0; i < nb ; i++) {
         for (j=0; j < na ; j++) b[j] = a[j][i+na];
         lubksb(a,na,indx,b);
         for (j=0; j < na ; j++) a[j][i+na] = b[j];
         }
 
      free(indx);
      free(b);
      }
