#include "file30.h"
#include "file30.gbl"

/*
** rd_max_am():  Reads in the maximum orbital quantum number of AOs in the basis.
**
**   takes no arguments.
**
**   returns: int max_am (0 corresponds to s-functions, 1 - to up to p-functions, etc.)
*/


int file30_rd_max_am(void)
{
  return info30_.mconst[8];
}
