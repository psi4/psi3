
/* $Log$
 * Revision 1.1  2000/02/04 22:53:21  evaleev
 * Initial revision
 *
/* Revision 2.1  1991/06/15 18:29:43  seidl
/* *** empty log message ***
/* */

static char *rcsid = "$Id$";

#include "includes.h"

void print_mat(a,m,n,out)
   double **a;
   int m,n;
   FILE *out;

   {
      int ii,jj,kk,nn,ll;
      int i,j,k;

      ii=0;jj=0;
L200:
      ii++;
      jj++;
      kk=10*jj;
      nn=n;
      if (nn > kk) nn=kk;
      ll = 2*(nn-ii+1)+1;
      fprintf (out,"\n");
      for (i=ii; i <= nn; i++) fprintf(out,"       %5d",i);
      fprintf (out,"\n");
      for (i=0; i < m; i++) {
         fprintf (out,"\n%5d",i+1);
         for (j=ii-1; j < nn; j++) {
            fprintf (out,"%12.7f",a[i][j]);
            }
         }
      fprintf (out,"\n");
      if (n <= kk) {
         fflush(out);
         return;
         }
      ii=kk; goto L200;
      }
