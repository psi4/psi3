/*!
** \file wreadw.c
*/

/* $Log$
 * Revision 1.2  2002/04/19 21:48:06  sherrill
 * Remove some unused functions and do doxygen markup of libciomr.
 *
/* Revision 1.1.1.1  2000/02/04 22:53:24  evaleev
/* Started PSI 3 repository
/*
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


/*
** wreadw: reads size bytes from tape into buffer starting at fword.
** nxtwrd is modified to give the new current pointer location after the
** read operation is completed.
**
**   \param tape   = file number
**   \param buffer = buffer to store read information
**   \param size   = number of bytes to be read
**   \param fword  = first byte of buffer to be read
**   \param nxtwrd = pointer to hold file pointer position after read
*/  
void wreadw(int tape, char *buffer, int size, PSI_FPTR fword, PSI_FPTR *nxtwrd)
   {
      iordr_(&tape,buffer,&fword,&size);
      *nxtwrd = fword + size;
      ptr.wptr[tape] = *nxtwrd;
   }
