
/* $Log$
 * Revision 1.1  2000/02/04 22:53:21  evaleev
 * Initial revision
 *
/* Revision 2.4  1999/11/01 20:10:58  evaleev
/* Added explicit extern declarations of functions within the library.
/*
/* Revision 2.3  1997/08/25 21:49:58  crawdad
/* Making changes for extension of PSI file size limit.
/*
 * Revision 2.2  1997/06/23  12:25:53  crawdad
 * Multiple changes to libciomr: Moved "param.h" to "iomrparam.h" to avoid
 *     conflicts with similarly named system file under linux.  Corrected type
 *    casting in rread(), rwrit(), sread(), and swrit() functions.  Corrected
 *     unclosed tmpdisks.dat file in sequential.c.  Corrected block_matrix() to
 *    avoid malloc'ing zero-length arrays.
 *
 * -Daniel
 *
 * Revision 2.1  1991/06/15  18:29:33  seidl
 * initial revision
 * */

static char *rcsid = "$Id$";

#include "iomrparam.h"
#include "includes.h"

extern int io_locate(FILE *, char[]);

int oldstyleinput(void)
{
  FILE *input;
  int ierr;

  input = fopen("input.dat","r");
  ierr = io_locate(input,"# FILES ##");
  fclose(input);
  if (ierr == 0) return(1);
  return 0;
  }
