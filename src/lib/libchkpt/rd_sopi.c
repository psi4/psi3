/*!
  \file rd_sopi.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_sopi():  Reads in the number of symmetry orbitals in each irrep.
**
**   takes no arguments.
**
**   returns:
**     int *sopi  an array which has an element for each irrep of the
**                 point group of the molecule (n.b. not just the ones
**                 with a non-zero number of basis functions). each 
**                 element contains the number of symmetry orbitals for
**                 that irrep. Also, see chkpt_rd_orbspi().
*/

int *chkpt_rd_sopi(void)
{
  int nirreps, *sopi;

  nirreps = chkpt_rd_nirreps();
  sopi = init_int_array(nirreps);
  psio_read_entry(PSIF_CHKPT, "::SO's per irrep", (char *) sopi, nirreps*sizeof(int));
  
  return sopi;
}
