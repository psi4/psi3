
/* $Id$ */
/* $Log$
 * Revision 1.1  2000/02/04 22:53:21  evaleev
 * Initial revision
 *
/* Revision 2.7  1997/09/12 13:52:54  crawdad
/* Changing marco name from ULL to PSI_FPTR.
/*
 * Revision 2.6  1997/08/25  21:49:59  crawdad
 * Making changes for extension of PSI file size limit.
 *
 * Revision 2.5  1995/04/01  20:52:58  fermann
 * changed bytewise file pointers such as first, last and length to long
 * unsigned ints in order to handle up to 4 gigabyte tmp files (striped into
 * individual pieces of less than 2 gigabytes).  added functions li2sec and
 * sec2li for where they are needed.
 *
 * Revision 2.4  1994/09/19  23:32:02  cljanss
 * Cleaned up allocation of globals.  Got rid of some globals and fixed a bug.
 *
 * Revision 2.3  1994/08/10  00:17:22  dcrawfrd
 * Added timing routines that are called by iomr and that give hostnames.
 *
 * Revision 2.2  1994/06/02  02:28:50  seidl
 * define ALLOC_GLOBALS
 *
 * Revision 2.1  1991/06/15  18:32:28  seidl
 * *** empty log message ***
 * */

#include "iomrparam.h"

#ifdef ALLOC_GLOBALS
#define EXTERN
#else
#define EXTERN extern
#endif

struct pointer {
     PSI_FPTR *wptr;
     };

EXTERN struct pointer ptr;
EXTERN int sector;
EXTERN time_t time_start, time_end;
#if !defined(SGI)
EXTERN struct tms total_tmstime;
#endif

#undef EXTERN
