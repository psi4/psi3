/*!
  \file add_arr.c
*/

/* $Log$
 * Revision 1.2  2002/03/25 02:43:45  sherrill
 * Update documentation
 *
/* Revision 1.1.1.1  2000/02/04 22:53:17  evaleev
/* Started PSI 3 repository
/*
/* Revision 2.1  1991/06/15 18:28:39  seidl
/* *** empty log message ***
/* */

static char *rcsid = "$Id$";

/*!
** add_arr: Add arrays a and b and put the result in array c.  Adds
** the first n elements
*/
void add_arr(double *a, double *b, double *c, int n)
   {
      register int i;

      for (i=0; i < n ; i++) {
         c[i] = a[i]+b[i];
         }
      }
