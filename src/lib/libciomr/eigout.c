/*!
  \file eigout.c
  \ingroup (CIOMR)
*/

/* $Log$
 * Revision 1.3  2002/06/01 18:23:54  sherrill
 * Upgrade doxygen documentation
 *
/* Revision 1.2  2002/04/18 21:47:35  sherrill
/* Here's some changes to document via doxygen and upgrade to ANSI C
/* instead of K&R declarations.
/*
/* Revision 1.1.1.1  2000/02/04 22:53:18  evaleev
/* Started PSI 3 repository
/*
/* Revision 2.3  1995/01/17 19:56:52  local
/* Fixed undeclared variable bug.  Minor modification.
/*
 * Revision 2.2  1995/01/16  22:48:32  cljanss
 * Minor changes to make the SGI compiler happy.
 *
 * Revision 2.1  1991/06/15  18:28:46  seidl
 * *** empty log message ***
 * */

static char *rcsid = "$Id$";

#include "includes.h"

/*!
** eigout(): Print out eigenvectors and eigenvalues.  Prints an n x m
** matrix of eigenvectors.  Under each eigenvector, the corresponding
** elements of two arrays, b and c, will also be printed.  This is
** useful for printing, for example, the SCF eigenvectors with their
** associated eigenvalues (orbital energies) and also the population.
**
** \ingroup (CIOMR)
*/
void eigout(double **a, double *b, double *c, int m, int n, FILE *out)
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
      fprintf (out,"\n     ");
      for (j=ii-1; j < nn; j++) {
         fprintf(out,"%12.7f",c[j]);
         }
      fprintf (out,"\n");
      if (n <= kk) {
         fflush(out);
         return;
         }
      ii=kk; goto L200;
      }
