/*!
   \file zero.c
   \ingroup (CIOMR)
*/

#include "includes.h"

/*!
** zero_arr: zero out an array of length 'size'.
** \ingroup (CIOMR)
*/
void zero_arr(double *a, int size)
{
  bzero(a,sizeof(double)*size);
}

/*!
** zero_mat: zero out a matrix 'a' with n rows and m columns 
** \ingroup (CIOMR)
*/
void zero_mat(double **a, int n, int m)
{
  register int i;

  for (i=0; i < n ; i++) {
    bzero(a[i],sizeof(double)*m);
  }
}

