/*!
  \file phase_check.c
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** int chkpt_rd_phase_check()
**
** Reads a boolean flag indicating whether the SCF code was able to correct
** the phases of the molecular orbitals relative to the guess orbitals.  This
** is important for restarting correlated wfn calculations from earlier vectors.
**
** arguments: none
**
** returns: int pcheck;
*/
int chkpt_rd_phase_check(void)
{
  int pcheck;

  psio_read_entry(PSIF_CHKPT, "::Phase check", (char *) &pcheck, sizeof(int));
  return pcheck;
}

/*!
** void chkpt_wt_phase_check(int)
**
** Reads a boolean flag indicating whether the SCF code was able to correct
** the phases of the molecular orbitals relative to the guess orbitals.  This
** is important for restarting correlated wfn calculations from earlier vectors.
**
** arguments: 
**  \param int pcheck
**
** returns: none
*/

void chkpt_wt_phase_check(int pcheck)
{
  psio_write_entry(PSIF_CHKPT, "::Phase check", (char *) &pcheck, sizeof(int));
}
