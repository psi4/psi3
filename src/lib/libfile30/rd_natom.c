#include "file30.h"
#include "file30.gbl"

/*
** file30_rd_natom():  Reads in the total number of atoms.
**
**   takes no arguments.
**
**   returns: int natom   total number of atoms.
*/

int file30_rd_natom(void)
{
  return info30_.mconst[18];
}
