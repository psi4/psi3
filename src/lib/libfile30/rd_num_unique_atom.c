/*!
  \file rd_num_unique_atom.c
*/

#include "file30.h"
#include "file30.gbl"

/*!
** file30_rd_num_unique_atom():  Reads in the number of symmetry unique atoms.
**
**   takes no arguments.
**
**   returns: int num_unique_atom   number of symmetry unique atoms.
*/

int file30_rd_num_unique_atom(void)
{
  return info30_.mconst[5];
}
