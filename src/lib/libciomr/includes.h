
/* $Id$ */
/* $Log$
 * Revision 1.1  2000/02/04 22:53:19  evaleev
 * Initial revision
 *
/* Revision 2.7  1997/06/23 12:25:47  crawdad
/* Multiple changes to libciomr: Moved "param.h" to "iomrparam.h" to avoid
/*     conflicts with similarly named system file under linux.  Corrected type
/*    casting in rread(), rwrit(), sread(), and swrit() functions.  Corrected
/*     unclosed tmpdisks.dat file in sequential.c.  Corrected block_matrix() to
/*    avoid malloc'ing zero-length arrays.
/*
/* -Daniel
/*
 * Revision 2.6  1995/03/09  00:07:27  psi
 * include sys/param.h on linux machines...otherwise user and sys times are
 * off by 100/60.
 *
 * Revision 2.5  1995/01/16  22:49:33  cljanss
 * Minor changes to make the SGI compiler happy.
 *
 * Revision 2.4  1994/08/10  00:17:16  dcrawfrd
 * Added timing routines that are called by iomr and that give hostnames.
 *
 * Revision 2.3  1994/06/02  02:29:35  seidl
 * add SGI specific includes
 *
 * Revision 2.2  1991/09/18  20:47:36  seidl
 * dec changes
 *
 * Revision 2.1  1991/06/15  18:32:26  seidl
 * *** empty log message ***
 * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include "libciomr.h"

#if defined(DEC)
#  include <sys/types.h>
#  include <time.h>
#  include <sys/times.h>
#  include <sys/param.h>
#  include <sys/mount.h>
#  include <sys/stat.h>
#elif defined(SGI)
#  include <sys/stat.h>
#  include <sys/types.h>
#  include <sys/times.h>
#  include <time.h>
#  include <sys/param.h>
#  include <sys/mount.h>
#elif defined(AIX)
#  include <sys/stat.h>
#  include <sys/vfs.h>
#  include <time.h>
#  include <sys/times.h>
#elif defined(Linux)
#  include <sys/stat.h>
#  include <sys/vfs.h>
#  include <time.h>
#  include <sys/times.h>
#  include <sys/param.h>
#else
#  include <sys/stat.h>
#  include <sys/vfs.h>
#  include <time.h>
#  include <sys/times.h>
#endif

#if (defined(DEC)||defined(SUN))
#define SEEK_SET 0
#define SEEK_CUR 1
#define SEEK_END 2
#endif

#define DSIGN(a,b) ((b) >= 0.0) ? (fabs(a)) : (-fabs(a))
#define MIN0(a,b) ((a)<(b)) ? (a) : (b)
#define MAX0(a,b) ((a)>(b)) ? (a) : (b)

#ifdef AIX
# define MALLOC malloc
#else
# define MALLOC malloc
#endif
