/*!
  \file rd_iopen.c
*/

#include "file30.h"
#include "file30.gbl"

/*!
** file30_rd_iopen():	Reads in dimensionality of ALPHA and BETA vectors 
**			of two-electron coupling coefficients for open shells.
**			Note : IOPEN = MM * (MM + 1), where MM is 
**			the total number of irreps containing singly 
**			occupied orbitals.
**
**  takes no arguments.
**
**  returns: int iopen	dimensionality of ALPHA and BETA vectors of coupling
**			coefficients for open shells. 
**			
*/


int file30_rd_iopen(void)
{
  return info30_.mconst[42];
}
