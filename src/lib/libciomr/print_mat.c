/*!
** \file print_mat.c
*/

/* $Log$
 * Revision 1.2  2002/04/19 21:48:06  sherrill
 * Remove some unused functions and do doxygen markup of libciomr.
 *
/* Revision 1.1.1.1  2000/02/04 22:53:21  evaleev
/* Started PSI 3 repository
/*
/* Revision 2.1  1991/06/15 18:29:43  seidl
/* *** empty log message ***
/* */

static char *rcsid = "$Id$";

#include "includes.h"

/*
** print_mat: Print a matrix a of dimensions mxn to file pointer out.
*/
void print_mat(double **a, int m, int n, FILE *out)
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
