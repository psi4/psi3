/*!
  \file orbspi.c
*/

#include <stdio.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_orbspi():  Reads in the number of molecular orbitals in each irrep.
**
**   takes no arguments.
**
**   returns:
**     int *orbspi  an array which has an element for each irrep of the
**                 point group of the molecule (n.b. not just the ones
**                 with a non-zero number of basis functions). each 
**                 element contains the number of molecular orbitals for
**                 that irrep. Also, see chkpt_rd_sopi().
*/

int *chkpt_rd_orbspi(void)
{
  int nirreps;
  int *orbspi;

  nirreps = chkpt_rd_nirreps();
  orbspi = init_int_array(nirreps);

  psio_read_entry(PSIF_CHKPT, "::MO's per irrep", (char *) orbspi, 
                  nirreps*sizeof(int));

  return orbspi;
}


/*!
** chkpt_wt_orbspi():  Writes the number of molecular orbitals in each irrep.
**
** \param orbspi = an array which has an element for each irrep of the
**                 point group of the molecule (n.b. not just the ones
**                 with a non-zero number of basis functions). each 
**                 element contains the number of molecular orbitals for
**                 that irrep. Also, see chkpt_rd_sopi().
**
** returns: none
*/

void chkpt_wt_orbspi(int *orbspi)
{
  int nirreps;

  nirreps = chkpt_rd_nirreps();

  psio_write_entry(PSIF_CHKPT, "::MO's per irrep", (char *) orbspi, 
                   nirreps*sizeof(int));
}
