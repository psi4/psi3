
/* $Log$
 * Revision 1.1  2000/02/04 22:53:17  evaleev
 * Initial revision
 *
/* Revision 2.1  1991/06/15 18:28:43  seidl
/* *** empty log message ***
/* */

static char *rcsid = "$Id$";

#include "includes.h"

void add_mat(a,b,c,n,m)
   double **a, **b, **c;
   int n,m;

   {
      register int i,j;

      if (n != m) {
         for (i=0; i < n ; i++) {
            for (j=0; j < m ; j++) {
               c[i][j] = a[i][j]+b[i][j];
               }
            }
         }
      else {
         for (i=0; i < n; i++) {
            for (j=0; j < i; j++) {
               c[i][j] = a[i][j]+b[i][j];
               c[j][i] = a[j][i]+b[j][i];
               }
            c[i][i] = a[i][i]+b[i][i];
            }
         }
      }
