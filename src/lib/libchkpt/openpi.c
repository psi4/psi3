/*!
  \file openpi.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_openpi():  Reads in the number of open-shell orbitals in each irrep.
**
**   takes no arguments.
**
**   returns:
**     int *openpi  an array which has an element for each irrep of the
**                 point group of the molecule (n.b. not just the ones
**                 with a non-zero number of basis functions). each 
**                 element contains the number of open-shell orbitals for
**                 that irrep.
*/

int *chkpt_rd_openpi(void)
{
  int nirreps;
  int *openpi;

  nirreps = chkpt_rd_nirreps();
  openpi = init_int_array(nirreps);

  psio_read_entry(PSIF_CHKPT, "::Open shells per irrep", (char *) openpi, nirreps*sizeof(int));

  return openpi;
}

/*!
** chkpt_wt_openpi():  Writes the number of open-shell orbitals in each irrep.
**
** arguments:
**  \param int *openpi  an array which has an element for each irrep of the
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

  nirreps = chkpt_rd_nirreps();

  psio_write_entry(PSIF_CHKPT, "::Open shells per irrep", (char *) openpi, nirreps*sizeof(int));
}

