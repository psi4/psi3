/*!
  \file ecorr.c
  \ingroup (CHKPT)
*/

#include <stdio.h>
#include <stdlib.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_ecorr():  Reads in the correlated energy.
**
** takes no arguments.
**
** returns: e_corr = the correlated energy.  To get some
**        information (a label) on the type of correlated wavefunction
**        used to get this energy, see rd_corr_lab().
** \ingroup (CHKPT)
*/
double chkpt_rd_ecorr(void)
{
  double ecorr;
  char *keyword;
  keyword = chkpt_build_keyword("Correlation energy");

  psio_read_entry(PSIF_CHKPT, keyword, (char *) &ecorr,
    sizeof(double));

  free(keyword);
  return ecorr;
}


/*!
** chkpt_wt_ecorr():  Writes out the correlated energy.
**
** \param e_corr = the correlated energy.  To get some
**        information (a label) on the type of correlated wavefunction
**        used to get this energy, see rd_corr_lab().
**
** returns: none
** \ingroup (CHKPT)
*/
void chkpt_wt_ecorr(double ecorr)
{
  char *keyword;
  keyword = chkpt_build_keyword("Correlation energy");

  psio_write_entry(PSIF_CHKPT, keyword, (char *) &ecorr,
    sizeof(double));

  free(keyword);
}

