
/* $Log$
 * Revision 1.2  2002/03/25 02:43:45  sherrill
 * Update documentation
 *
/* Revision 1.1.1.1  2000/02/04 22:53:17  evaleev
/* Started PSI 3 repository
/*
/* Revision 2.1  1991/06/15 18:28:43  seidl
/* *** empty log message ***
/* */

static char *rcsid = "$Id$";

#include "includes.h"

/*!
** add_mat(): Add matrices a and b into c for n rows and m columns
*/
void add_mat(double **a, double **b, double **c, int n, int m)
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
