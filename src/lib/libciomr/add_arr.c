
/* $Log$
 * Revision 1.1  2000/02/04 22:53:17  evaleev
 * Initial revision
 *
/* Revision 2.1  1991/06/15 18:28:39  seidl
/* *** empty log message ***
/* */

static char *rcsid = "$Id$";

void add_arr(a,b,c,n)
   double *a, *b, *c;
   int n;

   {
      register int i;

      for (i=0; i < n ; i++) {
         c[i] = a[i]+b[i];
         }
      }
