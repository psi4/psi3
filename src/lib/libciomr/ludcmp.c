
/* $Log$
 * Revision 1.2  2002/04/04 22:24:38  evaleev
 * Converted allocation functions (init_array, etc.) to take unsigned long ints
 * to be able to allocate properly 2GB+ chunks). Some declarations cleaned up.
 *
/* Revision 1.1.1.1  2000/02/04 22:53:20  evaleev
/* Started PSI 3 repository
/*
/* Revision 2.2  1999/11/01 20:10:57  evaleev
/* Added explicit extern declarations of functions within the library.
/*
/* Revision 2.1  1991/06/15 18:29:27  seidl
/* *** empty log message ***
/* */

static char *rcsid = "$Id$";

#include "includes.h"

void ludcmp(a,n,indx,d)
   int n, *indx;
   double **a, *d;

   {
      int i,imax,j,k;
      double big,dum,sum,temp;
      double *vv;

      vv = (double *) init_array(n);

      *d = 1.0;

      for (i=0; i < n ; i++) {
         big=0.0;
         for (j=0; j < n; j++) {
            if ((temp=fabs(a[i][j])) > big) big=temp;
            }
         if (big == 0.0) {
            *d = 0.0;
            return;
            }
         vv[i] = 1.0/big;
         }
      for (j=0; j < n ; j++) {
         for (i=0; i < j ; i++) {
            sum = a[i][j];
            for (k=0; k < i ; k++) sum -= a[i][k]*a[k][j];
            a[i][j] = sum;
            }
         big = 0.0;
         for (i=j ; i < n ; i++) {
            sum=a[i][j];
            for (k=0; k < j ; k++) sum -= a[i][k]*a[k][j];
            a[i][j] = sum;
            if ((dum=vv[i]*fabs(sum)) >= big) {
               big = dum;
               imax = i;
               }
            }
         if (j != imax) {
            for (k=0; k < n; k++) {
               dum=a[imax][k];
               a[imax][k]=a[j][k];
               a[j][k]=dum;
               }
            *d = -(*d);
            vv[imax]=vv[j];
            }
         indx[j]=imax;
         if (a[j][j] == 0.0) a[j][j] = 1.0e-20;
         if (j != n-1) {
            dum = 1.0/a[j][j];
            for (i=j+1; i < n ; i++) a[i][j] *= dum;
            }
         }
      free(vv);
      }
