/*!
  \file etot.c
*/

#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_etot(): Reads in the total energy.
**
**   takes no arguments.
**
**   returns: double etot  the total energy.
*/

double chkpt_rd_etot(void)
{
  double etot;
  psio_read_entry(PSIF_CHKPT, "::Total energy", (char *) &etot, sizeof(double));
  return etot;
}

/*!
** chkpt_wt_etot(): Writes out the total energy.
**
**  arguments: 
**   \param double etot  the total energy.
**
** returns: none
*/

void chkpt_wt_etot(double etot)
{
  psio_write_entry(PSIF_CHKPT, "::Total energy", (char *) &etot, sizeof(double));
}
