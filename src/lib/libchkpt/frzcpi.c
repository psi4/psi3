/*!
  \file frzcpi.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_frzcpi():  Reads in the number of frozen doubly occupied molecular orbitals in each irrep.
**
**   takes no arguments.
**
**   returns:
**     int *frzcpi  an array which has an element for each irrep of the
**                 point group of the molecule (n.b. not just the ones
**                 with a non-zero number of basis functions). each 
**                 element contains the number of frozen doubly occupied
**                 molecular orbitals for
**                 that irrep. Also, see chkpt_rd_sopi().
*/

int *chkpt_rd_frzcpi(void)
{
  int nirreps;
  int *frzcpi;
  char *keyword;
  keyword = chkpt_build_keyword("Frozen DOCC per irrep");

  nirreps = chkpt_rd_nirreps();
  frzcpi = init_int_array(nirreps);

  psio_read_entry(PSIF_CHKPT, keyword, (char *) frzcpi, nirreps*sizeof(int));

  free(keyword);
  return frzcpi;
}


/*!
** chkpt_wt_frzcpi():  Writes the number of frozen doubly occupied molecular orbitals in each irrep.
**
** \param frzcpi = an array which has an element for each irrep of the
**                 point group of the molecule (n.b. not just the ones
**                 with a non-zero number of basis functions). each 
**                 element contains the number of frozen doubly occupied molecular orbitals for
**                 that irrep. Also, see chkpt_rd_sopi().
**
** returns: none
*/

void chkpt_wt_frzcpi(int *frzcpi)
{
  int nirreps;
  char *keyword;
  keyword = chkpt_build_keyword("Frozen DOCC per irrep");

  nirreps = chkpt_rd_nirreps();

  psio_write_entry(PSIF_CHKPT, keyword, (char *) frzcpi, nirreps*sizeof(int));

  free(keyword);
}
