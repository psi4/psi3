/*!
  \file rd_disp.c
*/

#include "file30.h"
#include "file30.gbl"

/*!
** rd_disp():  Reads in the current geometry displacement number.
**
**   takes no arguments.
**
**   returns: int disp   the current geometry displacement number
*/


int file30_rd_disp(void)
{
  return info30_.mconst[51];
}
