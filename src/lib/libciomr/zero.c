/*!
   \file zero.c
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
/* Revision 2.1  1991/06/15 18:30:19  seidl
/* *** empty log message ***
/* */

static char *rcsid = "$Id$";

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

