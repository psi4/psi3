
/* $Log$
 * Revision 1.1  2000/02/04 22:53:19  evaleev
 * Initial revision
 *
/* Revision 2.4  1997/09/12 13:52:49  crawdad
/* Changing marco name from ULL to PSI_FPTR.
/*
 * Revision 2.3  1997/08/25  21:49:51  crawdad
 * Making changes for extension of PSI file size limit.
 *
 * Revision 2.2  1995/04/01  20:52:53  fermann
 * changed bytewise file pointers such as first, last and length to long
 * unsigned ints in order to handle up to 4 gigabyte tmp files (striped into
 * individual pieces of less than 2 gigabytes).  added functions li2sec and
 * sec2li for where they are needed.
 *
 * Revision 2.1  1991/06/15  18:29:16  seidl
 * *** empty log message ***
 * */

static char *rcsid = "$Id$";

/* converts integer file pointer to sector pointer */

#include "iomrparam.h"

int i2sec(n)
    PSI_FPTR n;
    {
      int num;

      num = (int) ((n+4096)/4096);
      return(num);
    }
