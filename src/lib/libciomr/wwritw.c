
/* $Log$
 * Revision 1.1  2000/02/04 22:53:24  evaleev
 * Initial revision
 *
/* Revision 2.5  1999/11/01 20:11:01  evaleev
/* Added explicit extern declarations of functions within the library.
/*
/* Revision 2.4  1997/09/12 13:53:09  crawdad
/* Changing marco name from ULL to PSI_FPTR.
/*
 * Revision 2.3  1997/08/25  21:50:14  crawdad
 * Making changes for extension of PSI file size limit.
 *
 * Revision 2.2  1995/04/01  20:53:12  fermann
 * changed bytewise file pointers such as first, last and length to long
 * unsigned ints in order to handle up to 4 gigabyte tmp files (striped into
 * individual pieces of less than 2 gigabytes).  added functions li2sec and
 * sec2li for where they are needed.
 *
 * Revision 2.1  1991/06/15  18:30:18  seidl
 * *** empty log message ***
 * */

static char *rcsid = "$Id$";

#include "iomrparam.h"
#include "includes.h"
#include "pointers.h"

extern void iowrr_(int *, char *, PSI_FPTR *, int *);

/* writes nwords bytes to unit from buffer starting at fword */
/* returns nxtwrd which is the current pointer location */

void wwritw(unit, buffer, nwords, fword, nxtwrd)
   int unit,nwords;
   char *buffer;
   PSI_FPTR fword, *nxtwrd;

   {
      iowrr_(&unit,buffer,&fword,&nwords);
      *nxtwrd = fword + nwords;
      ptr.wptr[unit]= *nxtwrd;
   }
