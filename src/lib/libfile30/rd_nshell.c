/*!
  \file rd_nshell.c
*/

#include "file30.h"
#include "file30.gbl"

/*!
** file30_rd_nshell():  Reads in the total number of shells. For example,
** DZP basis for carbon atom (9s/4s,5p/2p,1d/1d) has total 15 basis functions,
** 15 primitives, and 7 shells. Shells of all atoms are counted (compare nprim).
**
**   takes no arguments.
**
**   returns: int nshell   total number of shells.
*/


int file30_rd_nshell(void)
{
  return info30_.mconst[26];
}
