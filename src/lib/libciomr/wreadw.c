
/* $Log$
 * Revision 1.1  2000/02/04 22:53:24  evaleev
 * Initial revision
 *
/* Revision 2.5  1999/11/01 20:11:00  evaleev
/* Added explicit extern declarations of functions within the library.
/*
/* Revision 2.4  1997/09/12 13:53:08  crawdad
/* Changing marco name from ULL to PSI_FPTR.
/*
 * Revision 2.3  1997/08/25  21:50:13  crawdad
 * Making changes for extension of PSI file size limit.
 *
 * Revision 2.2  1995/04/01  20:53:11  fermann
 * changed bytewise file pointers such as first, last and length to long
 * unsigned ints in order to handle up to 4 gigabyte tmp files (striped into
 * individual pieces of less than 2 gigabytes).  added functions li2sec and
 * sec2li for where they are needed.
 *
 * Revision 2.1  1991/06/15  18:30:17  seidl
 * *** empty log message ***
 * */

static char *rcsid = "$Id$";

#include "iomrparam.h"
#include "includes.h"
#include "pointers.h"

extern void iordr_(int *, char *, PSI_FPTR *, int *);

/* reads size bytes from tape into buffer starting at fword */
/* returns nxtwrd, which is the current pointer location */

void wreadw(tape, buffer, size, fword, nxtwrd)
   int tape,size;
   PSI_FPTR fword, *nxtwrd;
   char *buffer;

   {
      iordr_(&tape,buffer,&fword,&size);
      *nxtwrd = fword + size;
      ptr.wptr[tape] = *nxtwrd;
   }
