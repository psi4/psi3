/*!
   \file tri_to_sq.c
   \ingroup (CIOMR)
*/

/* $Log$
 * Revision 1.3  2002/06/01 18:23:54  sherrill
 * Upgrade doxygen documentation
 *
/* Revision 1.2  2002/04/19 21:48:06  sherrill
/* Remove some unused functions and do doxygen markup of libciomr.
/*
/* Revision 1.1.1.1  2000/02/04 22:53:24  evaleev
/* Started PSI 3 repository
/*
/* Revision 2.1  1991/06/15 18:30:15  seidl
/* *** empty log message ***
/* */

static char *rcsid = "$Id$";

#include "includes.h"

/*!
** tri_to_sq: converts lower triangle to square matrix.
**
** \param amat = lower triangle matrix
** \param bmat = square matrix
** \param size = number of rows/cols of matrix
** \ingroup (CIOMR)
*/
void tri_to_sq(double *amat, double **bmat, int size)
{
  int i, j, ij;

  ij=0;
  for(i = 0 ; i < size ; i++) {
    for(j = 0 ; j <= i ; j++) {
      bmat[i][j] = amat[ij];
      bmat[j][i] = amat[ij++];
    }
  }
}

