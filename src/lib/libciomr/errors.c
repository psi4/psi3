
/* $Log$
 * Revision 1.1  2000/02/04 22:53:18  evaleev
 * Initial revision
 *
/* Revision 2.4  1999/11/01 20:10:55  evaleev
/* Added explicit extern declarations of functions within the library.
/*
/* Revision 2.3  1997/08/25 21:49:45  crawdad
/* Making changes for extension of PSI file size limit.
/*
 * Revision 2.2  1997/06/23  12:25:43  crawdad
 * Multiple changes to libciomr: Moved "param.h" to "iomrparam.h" to avoid
 *     conflicts with similarly named system file under linux.  Corrected type
 *    casting in rread(), rwrit(), sread(), and swrit() functions.  Corrected
 *     unclosed tmpdisks.dat file in sequential.c.  Corrected block_matrix() to
 *    avoid malloc'ing zero-length arrays.
 *
 * -Daniel
 *
 * Revision 2.1  1991/06/15  18:28:50  seidl
 * *** empty log message ***
 * */

static char *rcsid = "$Id$";


#include "iomrparam.h"
#include "includes.h"
#include "types.h"

void
no_path_given(name)
char *name;
{
  fprintf(stderr,"%s: no path given\n",name);
  ioabort();
  }

void
malloc_check(caller,data)
char *caller;
char *data;
{
  if (!data) {
    fprintf(stderr,"%s: malloc failed\n",caller);
    perror("malloc");
    ioabort();
    }
  }

void
fopen_check(caller,path,data)
char *caller;
char *path;
char *data;
{
  if (!data) {
    fprintf(stderr,"%s: fopen failed for %s\n",caller,path);
    perror("fopen");
    ioabort();
    }
  }

void
fread_error(caller)
char *caller;
{
  fprintf(stderr,"%s: fread failed\n",caller);
  perror("fread");
  ioabort();
  }

void
fwrite_error(caller)
char *caller;
{
  fprintf(stderr,"%s: fwrite failed\n",caller);
  perror("fwrite");
  ioabort();
  }
