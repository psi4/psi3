/*!
**\file sq_to_tri.c
*/

/* $Log$
 * Revision 1.2  2002/04/19 21:48:06  sherrill
 * Remove some unused functions and do doxygen markup of libciomr.
 *
/* Revision 1.1.1.1  2000/02/04 22:53:23  evaleev
/* Started PSI 3 repository
/*
/* Revision 2.1  1991/06/15 18:30:08  seidl
/* *** empty log message ***
/* */

static char *rcsid = "$Id$";

#include "includes.h"

/*
** sq_to_tri: converts square matrix to lower triangle
**
** \param bmat = matrix to convert
** \param amat = array to put lower triangle of bmat into
** \param size = number of rows/columns of bmat
*/ 
void sq_to_tri(double **bmat, double *amat, int size)
     {
        int i, j, ij;

        ij=0;
        for(i = 0 ; i < size ; i++) {
            for(j = 0 ; j <= i ; j++) {
               amat[ij++] = bmat[i][j];
               }
            }
     }

