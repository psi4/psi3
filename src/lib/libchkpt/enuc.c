/*!
  \file enuc.c
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_enuc(): Reads in the nuclear repulsion energy
**
**   takes no arguments.
**
**   returns: double enuc  the nuclear repulsion energy.
*/

double chkpt_rd_enuc(void)
{
  double enuc;

  psio_read_entry(PSIF_CHKPT, "::Nuclear rep. energy", (char *) &enuc, sizeof(double));
  return enuc;
}

/*!
** chkpt_wt_enuc(): Writes out the nuclear repulsion energy
**
** arguments: 
**  \param double enuc  the nuclear repulsion energy.
**
** returns: none
*/

void chkpt_wt_enuc(double enuc)
{
  psio_write_entry(PSIF_CHKPT, "::Nuclear rep. energy", (char *) &enuc, sizeof(double));
}
