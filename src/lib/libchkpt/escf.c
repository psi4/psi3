/*!
  \file escf.c
*/

#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_escf(): Reads in the scf energy.
**
**   takes no arguments.
**
**   returns: double escf  the scf energy.
*/

double chkpt_rd_escf(void)
{
  double escf;
  psio_read_entry(PSIF_CHKPT, "::SCF energy", (char *) &escf, sizeof(double));
  return escf;
}

/*!
** chkpt_wt_escf(): Writes out the scf energy.
**
**  arguments: 
**   \param double escf  the scf energy.
**
** returns: none
*/

void chkpt_wt_escf(double escf)
{
  psio_write_entry(PSIF_CHKPT, "::SCF energy", (char *) &escf, sizeof(double));
}
