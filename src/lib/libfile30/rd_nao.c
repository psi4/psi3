#include "file30.h"
#include "file30.gbl"

/*
** file30_rd_nao():  Reads in the total number of atomic orbitals.
**
**   takes no arguments.
**
**   returns: int nao   total number of atomic orbitals.
*/


int file30_rd_nao(void)
{
  return info30_.mconst[21];
}
