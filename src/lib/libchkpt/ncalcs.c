/*!
  \file calcs.c
  \ingroup (CHKPT)
*/

#include "chkpt.h"
#include <stdlib.h>
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
  char *keyword_mo, *keyword_alpha_mo;
  keyword_mo = chkpt_build_keyword("MO coefficients");
  keyword_alpha_mo = chkpt_build_keyword("Alpha MO coefficients");

  if (psio_tocscan(PSIF_CHKPT, keyword_mo) == NULL &&
      psio_tocscan(PSIF_CHKPT, keyword_alpha_mo) == NULL)
    return 0;
  else
    return 1;

  free(keyword_mo);
  free(keyword_alpha_mo);
}
