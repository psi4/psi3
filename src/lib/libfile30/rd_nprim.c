/*!
  \file rd_nprim.c
*/

#include "file30.h"
#include "file30.gbl"

/*!
** file30_rd_nprim():  Reads in the total number of primitive Gaussian 
**                     functions (only primitives of symmetry independent 
**                     atoms are taken into account!).
**
**   takes no arguments.
**
**   returns: int nprim   total number of primitive Gaussian functions.
*/


int file30_rd_nprim(void)
{
  return info30_.mconst[31];
}
