/*!
  \file rd_nmo.c
*/

#include "file30.h"
#include "file30.gbl"

/*!
** rd_nmo():  Reads in the total number of molecular orbitals.
**
**   takes no arguments.
**
**   returns: int nmo   total number of molecular orbitals.
*/


int file30_rd_nmo(void)
{
  return info30_.mconst[45];
}
