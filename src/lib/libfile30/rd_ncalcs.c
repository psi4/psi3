#include "file30.h"
#include "file30.gbl"

/*
** file30_rd_ncalcs():  Reads in the total number of calculations.
**
**   takes no arguments.
**
**   returns: int ncalcs   total number of calculations in file30.
*/

int file30_rd_ncalcs(void)
{
  return info30_.mconst[44];
}
