#include "file30.h"
#include "file30.gbl"

/*
** rd_nso():  Reads in the total number of SOs.
**
**   takes no arguments.
**
**   returns: int nso   total number of symmetry-adapted basis functions.
*/


int file30_rd_nso(void)
{
  return info30_.mconst[17];
}
