/*!
  \file enuc.c
  \ingroup (CHKPT)
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
**
** \ingroup (CHKPT)
*/

double chkpt_rd_enuc(void)
{
  double enuc;
  char *key;

  key = chkpt_build_keyword("Nuclear rep. energy");

  psio_read_entry(PSIF_CHKPT, key, (char *) &enuc, sizeof(double));
  free(key);
  return enuc;
}

/*!
** chkpt_wt_enuc(): Writes out the nuclear repulsion energy
**
** \param enuc = the nuclear repulsion energy.
**
** returns: none
** \ingroup (CHKPT)
*/

void chkpt_wt_enuc(double enuc)
{
  char *key;
  key = chkpt_build_keyword("Nuclear rep. energy");

  psio_write_entry(PSIF_CHKPT, key, (char *) &enuc, sizeof(double));
  free(key);
}
