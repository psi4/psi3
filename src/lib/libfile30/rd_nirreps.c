#include "file30.h"
#include "file30.gbl"

/*
** rd_nirreps():  Reads in the total number of irreducible representations
**   in the point group in which the molecule is being considered.
**
**   takes no arguments.
**
**   returns: int nirreps   total number of irreducible representations.
*/

int file30_rd_nirreps(void)
{
  return info30_.mconst[27];
}
