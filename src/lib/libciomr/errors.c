/*!
  \file errors.c
  \brief Print error messages and abort for various errors
  \ingroup (CIOMR)
*/

/* $Log$
 * Revision 1.3  2002/06/01 18:23:54  sherrill
 * Upgrade doxygen documentation
 *
/* Revision 1.2  2002/04/18 21:47:35  sherrill
/* Here's some changes to document via doxygen and upgrade to ANSI C
/* instead of K&R declarations.
/*
/* Revision 1.1.1.1  2000/02/04 22:53:18  evaleev
/* Started PSI 3 repository
/*
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

/*!
** no_path_given(): Print error message for no path given and abort
**
** \param name = name of calling routine
**
** \ingroup (CIOMR)
*/ 
void no_path_given(char *name)
{
  fprintf(stderr,"%s: no path given\n",name);
  ioabort();
}

/*!
** malloc_check(): Check to see if malloc succeeded or failed.  If failure,
** print error and abort.
**
** \param caller = name of calling routine
** \param data = pointer to new data (supposed to be char *, not very
**               useful anymore...) 
** \ingroup (CIOMR)
*/
void malloc_check(char *caller, char *data)
{
  if (!data) {
    fprintf(stderr,"%s: malloc failed\n",caller);
    perror("malloc");
    ioabort();
    }
}

/*!
** fopen_check(): See if fopen worked; if not, print error and abort
** \param caller = name of calling routine
** \param path = path for fopen
** \param data = pointer for output stream (probably shouldn't really 
**               be char *)
** \ingroup (CIOMR)
*/
void fopen_check(char *caller, char *path, char *data)
{
  if (!data) {
    fprintf(stderr,"%s: fopen failed for %s\n",caller,path);
    perror("fopen");
    ioabort();
    }
}

/*!
** fread_error(): If error in fread, print error and abort
** \ingroup (CIOMR)
*/
void fread_error(char *caller)
{
  fprintf(stderr,"%s: fread failed\n",caller);
  perror("fread");
  ioabort();
}

/*!
** fwrite_error(): If error in fwrite, print error and abort
** \ingroup (CIOMR)
*/
void fwrite_error(char *caller)
{
  fprintf(stderr,"%s: fwrite failed\n",caller);
  perror("fwrite");
  ioabort();
}
