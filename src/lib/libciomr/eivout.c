
/* $Log$
 * Revision 1.1  2000/02/04 22:53:18  evaleev
 * Initial revision
 *
/* Revision 2.2  1995/01/16 22:48:46  cljanss
/* Minor changes to make the SGI compiler happy.
/*
 * Revision 2.1  1991/06/15  18:28:48  seidl
 * *** empty log message ***
 * */

static char *rcsid = "$Id$";

#include "includes.h"

void eivout(a,b,m,n,out)
   double **a, *b;
   int m,n;
   FILE *out;

   {
      int ii,jj,kk,nn;
      int i,j;

      ii=0;jj=0;
L200:
      ii++;
      jj++;
      kk=10*jj;
      nn=n;
      if (nn > kk) nn=kk;
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
      fprintf (out,"\n     ");
      for (j=ii-1; j < nn; j++) {
         fprintf(out,"%12.7f",b[j]);
         }
      fprintf (out,"\n");
      if (n <= kk) {
         fflush(out);
         return;
         }
      ii=kk; goto L200;
      }
