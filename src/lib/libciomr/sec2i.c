
/* $Log$
 * Revision 1.1  2000/02/04 22:53:22  evaleev
 * Initial revision
 *
/* Revision 2.4  1997/09/12 13:53:02  crawdad
/* Changing marco name from ULL to PSI_FPTR.
/*
 * Revision 2.3  1997/08/25  21:50:07  crawdad
 * Making changes for extension of PSI file size limit.
 *
 * Revision 2.2  1995/04/01  20:53:04  fermann
 * changed bytewise file pointers such as first, last and length to long
 * unsigned ints in order to handle up to 4 gigabyte tmp files (striped into
 * individual pieces of less than 2 gigabytes).  added functions li2sec and
 * sec2li for where they are needed.
 *
 * Revision 2.1  1991/06/15  18:30:02  seidl
 * *** empty log message ***
 * */

static char *rcsid = "$Id$";

#include "includes.h"
#include "pointers.h"

/* converts sector pointer to integer pointer */
/* returns address of 1st element of sector n */

PSI_FPTR sec2i(n)
   int n;
   {
      PSI_FPTR num;
  
      num=(PSI_FPTR) (sizeof(int)*sector*n);
      return(num);
    }
