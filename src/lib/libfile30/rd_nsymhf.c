/*!
  \file rd_nsymhf.c
*/

#include "file30.h"
#include "file30.gbl"

/*!
** file30_rd_nsymhf():  Reads in the total number of irreps
**   in the point group in which the molecule is being considered which 
**   have non-zero number of basis functions.
**
**   takes no arguments.
**
**   returns: int nirreps   total number of irreducible representations
**      with a non-zero number of basis functions. For STO or DZ water, for
**      example, this is three, even though nirreps is 4 (see rd_nirreps()).
*/

int file30_rd_nsymhf(void)
{
  return info30_.mconst[40];
}
