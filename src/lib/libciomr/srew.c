
/* $Log$
 * Revision 1.1  2000/02/04 22:53:23  evaleev
 * Initial revision
 *
/* Revision 2.1  1991/06/15 18:30:10  seidl
/* *** empty log message ***
/* */

static char *rcsid = "$Id$";

#include "includes.h"
#include "pointers.h"

void srew(tape)
   int tape;

   {
      ptr.wptr[tape] = 0;
      }
