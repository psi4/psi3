
/* $Log$
 * Revision 1.2  2002/06/01 18:23:54  sherrill
 * Upgrade doxygen documentation
 *
/* Revision 1.1.1.1  2000/02/04 22:53:24  evaleev
/* Started PSI 3 repository
/*
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

/*!
** \file wwritw.c
** \ingroup (CIOMR)
*/

#include "iomrparam.h"
#include "includes.h"
#include "pointers.h"

extern void iowrr_(int *, char *, PSI_FPTR *, int *);

/*!
** wwritw: writes nwords bytes to unit from buffer starting at fword
** and returns nxtwrd which is the current pointer location.
**
** \param unit   = file number
** \param buffer = buffer holding info to write
** \param nwords = number of bytes to write
** \param fword  = file pointer to where to put first byte written
** \param nxtwrd = file pointer to next byte on disk (returned value)
**
** \ingroup (CIMOR)
*/
void wwritw(int unit, char *buffer, int nwords, PSI_FPTR fword, 
            PSI_FPTR *nxtwrd)
{
  iowrr_(&unit,buffer,&fword,&nwords);
  *nxtwrd = fword + nwords;
  ptr.wptr[unit]= *nxtwrd;
}

