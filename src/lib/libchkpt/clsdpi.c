/*!
  \file clsdpi.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_clsdpi():  Reads in the number of closed-shell orbitals in each irrep.
**
**   takes no arguments.
**
**   returns:
**     int *clsdpi  an array which has an element for each irrep of the
**                 point group of the molecule (n.b. not just the ones
**                 with a non-zero number of basis functions). each 
**                 element contains the number of closed-shell orbitals for
**                 that irrep.
*/

int *chkpt_rd_clsdpi(void)
{
  int nirreps;
  int *clsdpi;

  nirreps = chkpt_rd_nirreps();
  clsdpi = init_int_array(nirreps);

  psio_read_entry(PSIF_CHKPT, "::Closed shells per irrep", (char *) clsdpi, nirreps*sizeof(int));

  return clsdpi;
}


/*!
** chkpt_wt_clsdpi():  Writes the number of closed-shell orbitals in each irrep.
**
**  arguments:
**   \param int *clsdpi  an array which has an element for each irrep of the
**                 point group of the molecule (n.b. not just the ones
**                 with a non-zero number of basis functions). each 
**                 element contains the number of closed-shell orbitals for
**                 that irrep.
**
** returns: none
*/

void chkpt_wt_clsdpi(int *clsdpi)
{
  int nirreps;

  nirreps = chkpt_rd_nirreps();

  psio_write_entry(PSIF_CHKPT, "::Closed shells per irrep", (char *) clsdpi, nirreps*sizeof(int));
}
