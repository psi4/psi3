/*!
  \file frzvpi.c
*/

#include <stdio.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_frzvpi():  Reads in the number of frozen unoccupied molecular orbitals in each irrep.
**
**   takes no arguments.
**
**   returns:
**     int *frzvpi  an array which has an element for each irrep of the
**                 point group of the molecule (n.b. not just the ones
**                 with a non-zero number of basis functions). each 
**                 element contains the number of frozen unoccupied
**                 molecular orbitals for
**                 that irrep. Also, see chkpt_rd_sopi().
*/

int *chkpt_rd_frzvpi(void)
{
  int nirreps;
  int *frzvpi;

  nirreps = chkpt_rd_nirreps();
  frzvpi = init_int_array(nirreps);

  psio_read_entry(PSIF_CHKPT, "::Frozen UOCC per irrep", (char *) frzvpi, 
                  nirreps*sizeof(int));

  return frzvpi;
}


/*!
** chkpt_wt_frzvpi():  Writes the number of frozen unoccupied molecular orbitals in each irrep.
**
** \param frzvpi = an array which has an element for each irrep of the
**                 point group of the molecule (n.b. not just the ones
**                 with a non-zero number of basis functions). each 
**                 element contains the number of frozen unoccupied molecular orbitals for
**                 that irrep. Also, see chkpt_rd_sopi().
**
** returns: none
*/

void chkpt_wt_frzvpi(int *frzvpi)
{
  int nirreps;

  nirreps = chkpt_rd_nirreps();

  psio_write_entry(PSIF_CHKPT, "::Frozen UOCC per irrep", (char *) frzvpi, 
                   nirreps*sizeof(int));
}
