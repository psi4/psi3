/*!
  \file eivout.c
  \ingroup (CIOMR)
*/
 
/* $Log$
 * Revision 1.3  2002/06/01 18:23:54  sherrill
 * Upgrade doxygen documentation
 *
/* Revision 1.2  2002/04/19 21:48:06  sherrill
/* Remove some unused functions and do doxygen markup of libciomr.
/*
/* Revision 1.1.1.1  2000/02/04 22:53:18  evaleev
/* Started PSI 3 repository
/*
/* Revision 2.2  1995/01/16 22:48:46  cljanss
/* Minor changes to make the SGI compiler happy.
/*
 * Revision 2.1  1991/06/15  18:28:48  seidl
 * *** empty log message ***
 * */

static char *rcsid = "$Id$";

#include "includes.h"

/*!
** eivout: Print out eigenvectors and eigenvalues to the output file
**
** \param a = eigenvectors
** \param b = eigenvalues
** \param m = rows of a
** \param n = columns of a
** \param out = output file pointer
**
** \ingroup (CIOMR)
*/
void eivout(double **a, double *b, int m, int n, FILE *out)
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
