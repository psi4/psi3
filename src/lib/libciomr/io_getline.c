
/* $Log$
 * Revision 1.1  2000/02/04 22:53:19  evaleev
 * Initial revision
 *
/* Revision 2.2  1997/06/23 12:25:49  crawdad
/* Multiple changes to libciomr: Moved "param.h" to "iomrparam.h" to avoid
/*     conflicts with similarly named system file under linux.  Corrected type
/*    casting in rread(), rwrit(), sread(), and swrit() functions.  Corrected
/*     unclosed tmpdisks.dat file in sequential.c.  Corrected block_matrix() to
/*    avoid malloc'ing zero-length arrays.
/*
/* -Daniel
/*
 * Revision 2.1  1991/06/15  18:29:23  seidl
 * *** empty log message ***
 * */

static char *rcsid = "$Id$";


#include <stdio.h>
#include "iomrparam.h"

int
io_getline(input,line)
FILE *input;
char line[MAX_STRING];
{
  int i;

  for (i=0; i<MAX_STRING; i++) {
    if ((line[i]=fgetc(input)) == EOF) return(-1);
    if(feof(input)) return(-1);
    if (line[i] == '\n') {
      line[i] = '\0';
      return(0);
      }
    }

  fprintf(stdout,"io_getline: buffer size exceeded\n");
  fprintf(stderr,"io_getline: buffer size exceeded\n");
  return(-1);
  }
