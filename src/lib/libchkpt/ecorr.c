/*!
  \file ecorr.c
*/

#include <stdio.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_ecorr():  Reads in the correlated energy.
**
**    takes no arguments.
**
**    returns: double e_corr   the correlated energy.  To get some
**        information (a label) on the type of correlated wavefunction
**        used to get this energy, see rd_corr_lab().
*/
double chkpt_rd_ecorr(void)
{
  double ecorr;

  psio_read_entry(PSIF_CHKPT, "::Correlation energy", (char *) &ecorr,
		  sizeof(double));

  return ecorr;
}

/*!
** chkpt_wt_ecorr():  Writes out the correlated energy.
**
**  arguments: 
**   \param double e_corr   the correlated energy.  To get some
**        information (a label) on the type of correlated wavefunction
**        used to get this energy, see rd_corr_lab().
**
** returns: none
*/
void chkpt_wt_ecorr(double ecorr)
{
  psio_write_entry(PSIF_CHKPT, "::Correlation energy", (char *) &ecorr,
		   sizeof(double));
}

