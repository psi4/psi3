/*!
  \file phase_check.c
  \ingroup (CHKPT)
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
** returns: pcheck = Phase check flag (1 if phase has been checked, else 0)
**
** \ingroup (CHKPT)
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
** \param pcheck = Phase check flag (1 if phase has been checked, else 0)
**
** returns: none
** 
** \ingroup (CHKPT)
*/

void chkpt_wt_phase_check(int pcheck)
{
  psio_write_entry(PSIF_CHKPT, "::Phase check", (char *) &pcheck, sizeof(int));
}
