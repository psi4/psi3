
/* $Log$
 * Revision 1.1  2000/02/04 22:53:18  evaleev
 * Initial revision
 *
/* Revision 2.2  1998/02/03 19:34:07  evaleev
/* Modified eigsort(), rsp(), and sq_rsp() to sort eigenvalues and
/* eigenvectors in either ascending OR descending order.
/*
 * Revision 2.1  1991/06/15  18:28:47  seidl
 * *** empty log message ***
 * */

static char *rcsid = "$Id$";

#include "includes.h"

void eigsort(d,v,n)
   double *d, **v;
   int n;

   {
      int i,j,k;
      double p;

/* Modified by Ed - if n is negative - sort eigenvalues in descending order */

      if (n >= 0) {
        for (i=0; i < n-1 ; i++) {
           k=i;
           p=d[i];
           for (j=i+1; j < n; j++) {
              if (d[j] < p) {
                 k=j;
                 p=d[j];
                 }
              }
           if (k != i) {
              d[k]=d[i];
              d[i]=p;
              for (j=0; j < n; j++) {
                 p=v[j][i];
                 v[j][i]=v[j][k];
                 v[j][k]=p;
                 }
               }
            }
        }
      else {
        n = abs(n);
        for (i=0; i < n-1 ; i++) {
           k=i;
           p=d[i];
           for (j=i+1; j < n; j++) {
              if (d[j] > p) {
                 k=j;
                 p=d[j];
                 }
              }
           if (k != i) {
              d[k]=d[i];
              d[i]=p;
              for (j=0; j < n; j++) {
                  p=v[j][i];
                  v[j][i]=v[j][k];
                  v[j][k]=p;
                  }
              }
           }
        }
   }
