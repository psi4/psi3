/*!
  \file mxcoef.c
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** int chkpt_rd_mxcoef()  
** Reads the value of the constant mxcoef.
**  
**  returns: int mxcoef   the sum of the squares of the number of symmetry
**    orbitals for each irrep.  This gives the number of elements in the
**    non-zero symmetry blocks of the SCF eigenvector.  For STO water 
**    mxcoef = (4*4) + (0*0) + (1*1) + (2*2) = 21.  
*/

int chkpt_rd_mxcoef(void)
{
  int mxcoef;

  psio_read_entry(PSIF_CHKPT, "::Mxcoef", (char *) &mxcoef, sizeof(int));
  return mxcoef;
}

/*!
** void chkpt_wt_mxcoef(int)  
** Writes out the value of the constant mxcoef.
**  
**  arguments: 
**   \param int mxcoef   the sum of the squares of the number of symmetry
**    orbitals for each irrep.  This gives the number of elements in the
**    non-zero symmetry blocks of the SCF eigenvector.  For STO water 
**    mxcoef = (4*4) + (0*0) + (1*1) + (2*2) = 21.  
*/

void chkpt_wt_mxcoef(int mxcoef)
{
  psio_write_entry(PSIF_CHKPT, "::Mxcoef", (char *) &mxcoef, sizeof(int));
}

