#include "file30.h"
#include "file30.gbl"

/*
** file30_rd_nentry():  If zmatrix, reads total number of entries (including dummy atoms)
**
**   takes no arguments.
**
**   returns: int nentry total number of entries.
*/

int file30_rd_nentry(void)
{
  return info30_.mconst[19];
}
