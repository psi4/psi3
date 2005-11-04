/*!
  \file cdsalc2cd.c
  \ingroup (CHKPT)
*/

#include <stdio.h>
#include <stdlib.h>
#include <psifiles.h>
#include "chkpt.h"
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_cdsalcpi(): Read in number of SALCs per irrep
**
** takes no arguments.
**
** returns: cdsalcpi = An array of nirreps integers.
** 
** \ingroup (CHKPT)
*/

int *chkpt_rd_cdsalcpi(void)
{
  const int nirreps = chkpt_rd_nirreps();
  int *cdsalcpi = init_int_array(nirreps);
  psio_address ptr = PSIO_ZERO;
  char *keyword = chkpt_build_keyword("cartdisp SALCs per irrep");

  psio_read(PSIF_CHKPT, keyword, (char *) cdsalcpi, nirreps*sizeof(int), ptr, &ptr);

  free(keyword);
  return cdsalcpi;
}


/*!
** chkpt_wt_cdsalcpi(): Writes out number of SALCs per irrep
**
** \param cdsalcpi = An array of nirreps integers
**
** returns: none
**
** \ingroup (CHKPT)
*/

void chkpt_wt_cdsalcpi(const int *cdsalcpi)
{
  const int nirreps = chkpt_rd_nirreps();
  psio_address ptr = PSIO_ZERO;
  char *keyword = chkpt_build_keyword("cartdisp SALCs per irrep");

  psio_write(PSIF_CHKPT, keyword, (char *) cdsalcpi, nirreps*sizeof(int), ptr, &ptr);

  free(keyword);
}

