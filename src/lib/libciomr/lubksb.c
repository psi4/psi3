#include "includes.h"

void lubksb(a,n,indx,b)
   double **a, *b;
   int n, *indx;

   {
      int i,ii,ip,j;
      int t=0;
      double sum;

      for (i=0; i < n ; i++) {
         ip = indx[i];
         sum = b[ip];
         b[ip]=b[i];
         if(t) {
            for (j=ii; j <= i-1 ; j++) sum -= a[i][j]*b[j];
            }
         else if(sum) {
            ii=i;
            t++;
            }
         b[i]=sum;
         }
      for (i=n-1; i >= 0 ; i--) {
         sum = b[i];
         for (j=i+1; j < n ; j++) sum -= a[i][j]*b[j];
         b[i] = sum/a[i][i];
         }
      }
