/*!
  \file openpi.c
  \ingroup (CHKPT)
*/

#include <stdio.h>
#include <stdlib.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_openpi(): Reads in the number of open-shell orbitals in each irrep.
**
**   takes no arguments.
**
**   returns:
**     *openpi  an array which has an element for each irrep of the
**              point group of the molecule (n.b. not just the ones
**              with a non-zero number of basis functions). each 
**              element contains the number of open-shell orbitals for
**              that irrep.
*/

int *chkpt_rd_openpi(void)
{
  int nirreps;
  int *openpi;
  char *keyword;
  keyword = chkpt_build_keyword("Open shells per irrep");

  nirreps = chkpt_rd_nirreps();
  openpi = init_int_array(nirreps);

  psio_read_entry(PSIF_CHKPT, keyword, (char *) openpi, 
    nirreps*sizeof(int));

  free(keyword);
  return openpi;
}


/*!
** chkpt_wt_openpi():  Writes the number of open-shell orbitals in each irrep.
**
** \param *openpi = an array which has an element for each irrep of the
**                 point group of the molecule (n.b. not just the ones
**                 with a non-zero number of basis functions). each 
**                 element contains the number of open-shell orbitals for
**                 that irrep.
**
** returns: none
*/

void chkpt_wt_openpi(int *openpi)
{
  int nirreps;
  char *keyword;
  keyword = chkpt_build_keyword("Open shells per irrep");

  nirreps = chkpt_rd_nirreps();

  psio_write_entry(PSIF_CHKPT, keyword, (char *) openpi, 
    nirreps*sizeof(int));

  free(keyword);
}

