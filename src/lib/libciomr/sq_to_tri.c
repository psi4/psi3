/*!
   \file sq_to_tri.c
   \ingroup (CIOMR)
*/

#include "includes.h"

/*
** sq_to_tri: converts square matrix to lower triangle
**
** \param bmat = matrix to convert
** \param amat = array to put lower triangle of bmat into
** \param size = number of rows/columns of bmat
**
** \ingroup (CIOMR)
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

