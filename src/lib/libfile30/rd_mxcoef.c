#include "file30.h"
#include "file30.gbl"

/* file30_rd_mxcoef:  Reads the value of the constant mxcoef.
**  
**  takes no arguments.
** 
**  returns: int mxcoef   the sum of the squares of the number of symmetry
**    orbitals for each irrep.  This gives the number of elements in the
**    non-zero symmetry blocks of the SCF eigenvector.  For STO water 
**    mxcoef = (4*4) + (0*0) + (1*1) + (2*2) = 21.  
*/

int file30_rd_mxcoef(void)
{
  return info30_.mconst[41];
}
