
/* $Log$
 * Revision 1.1  2000/02/04 22:53:22  evaleev
 * Initial revision
 *
/* Revision 2.4  1997/09/12 13:52:57  crawdad
/* Changing marco name from ULL to PSI_FPTR.
/*
 * Revision 2.3  1997/08/25  21:50:02  crawdad
 * Making changes for extension of PSI file size limit.
 *
 * Revision 2.2  1995/04/01  20:52:59  fermann
 * changed bytewise file pointers such as first, last and length to long
 * unsigned ints in order to handle up to 4 gigabyte tmp files (striped into
 * individual pieces of less than 2 gigabytes).  added functions li2sec and
 * sec2li for where they are needed.
 *
 * Revision 2.1  1991/06/15  18:29:53  seidl
 * *** empty log message ***
 * */

static char *rcsid = "$Id$";

#include "iomrparam.h"
#include "includes.h"
#include "pointers.h"

void rgetsa(unit,iadr)
   int unit;
   int *iadr;

   {
      PSI_FPTR ipos,irec,test;
    
      ipos=ptr.wptr[unit];
      irec=(ipos/sizeof(PSI_FPTR)*sector)+1;
      test=sizeof(PSI_FPTR)*sector*(irec-1);

      if (ipos != test) {
         fprintf(stderr,"error encountered in rgetsa for file %d",unit);
         fprintf(stderr,"ipos,test = %lu %lu",ipos,test);
         exit(unit);
         }

      *iadr=irec;
  }
