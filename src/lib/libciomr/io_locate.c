
/* $Log$
 * Revision 1.1  2000/02/04 22:53:19  evaleev
 * Initial revision
 *
/* Revision 2.3  1999/11/01 20:10:56  evaleev
/* Added explicit extern declarations of functions within the library.
/*
/* Revision 2.2  1997/06/23 12:25:50  crawdad
/* Multiple changes to libciomr: Moved "param.h" to "iomrparam.h" to avoid
/*     conflicts with similarly named system file under linux.  Corrected type
/*    casting in rread(), rwrit(), sread(), and swrit() functions.  Corrected
/*     unclosed tmpdisks.dat file in sequential.c.  Corrected block_matrix() to
/*    avoid malloc'ing zero-length arrays.
/*
/* -Daniel
/*
 * Revision 2.1  1991/06/15  18:29:24  seidl
 * *** empty log message ***
 * */

static char *rcsid = "$Id$";


#include "includes.h"
#include "iomrparam.h"

extern int io_getline(FILE *, char[]);

int
io_locate(input,loc_token)
FILE *input;
char loc_token[MAX_STRING];
{
  char line[MAX_STRING];

  fseek(input,0L,0);

  for (;;) {
    if (io_getline(input,line) != 0) return(-1);
    if (!strncmp(loc_token,line,10)) {
      return(0);
      }
    }
  }

