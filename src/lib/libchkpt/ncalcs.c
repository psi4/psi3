/*!
  \file calcs.c
  \ingroup (CHKPT)
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** int chkpt_rd_ncalcs()  
** Reads in the total number of HF wave functions.
**
** returns: ncalcs = total number of HF wave functions in checkpoint
** \ingroup (CHKPT)
*/

int chkpt_rd_ncalcs(void)
{
  if (psio_tocscan(PSIF_CHKPT,"::MO coefficients") == NULL &&
      psio_tocscan(PSIF_CHKPT,"::Alpha MO coefficients") == NULL)
    return 0;
  else
    return 1;
}
