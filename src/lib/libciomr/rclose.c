/*!
** \file rclose.c
*/

/* $Log$
 * Revision 1.2  2002/04/19 21:48:06  sherrill
 * Remove some unused functions and do doxygen markup of libciomr.
 *
/* Revision 1.1.1.1  2000/02/04 22:53:22  evaleev
/* Started PSI 3 repository
/*
/* Revision 2.3  1999/11/01 20:10:58  evaleev
/* Added explicit extern declarations of functions within the library.
/*
/* Revision 2.2  1997/06/23 12:25:56  crawdad
/* Multiple changes to libciomr: Moved "param.h" to "iomrparam.h" to avoid
/*     conflicts with similarly named system file under linux.  Corrected type
/*    casting in rread(), rwrit(), sread(), and swrit() functions.  Corrected
/*     unclosed tmpdisks.dat file in sequential.c.  Corrected block_matrix() to
/*    avoid malloc'ing zero-length arrays.
/*
/* -Daniel
/*
 * Revision 2.1  1991/06/15  18:29:47  seidl
 * *** empty log message ***
 * */

static char *rcsid = "$Id$";

extern void ioclos_(int *, int *);

/*!
** rclose: close a binary file
**
** \param unit = file number
** \param status = 3 to keep file, 4 to erase it
*/
void rclose(int unit, int status)
{
      ioclos_(&unit,&status);
   }
