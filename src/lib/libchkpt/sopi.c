/*!
  \file sopi.c
  \ingroup (CHKPT)
*/

#include <stdio.h>
#include <stdlib.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_sopi()
** Reads in the number of symmetry orbitals in each irrep.
**
**  takes no arguments.
**
**  returns:
**    sopi =  an array which has an element for each irrep of the
**            point group of the molecule (n.b. not just the ones
**            with a non-zero number of basis functions). each 
**            element contains the number of symmetry orbitals for
**            that irrep. Also, see chkpt_rd_orbspi().
**
** \ingroup (CHKPT)
*/

int *chkpt_rd_sopi(void)
{
  int nirreps, *sopi;
  char *keyword;
  keyword = chkpt_build_keyword("SO's per irrep");

  nirreps = chkpt_rd_nirreps();
  sopi = init_int_array(nirreps);

  psio_read_entry(PSIF_CHKPT, keyword, (char *) sopi, nirreps*sizeof(int));
  
  free(keyword);
  return sopi;
}


/*!
** chkpt_wt_sopi():  Writes out the number of symmetry orbitals in each irrep.
**
** \param sopi = an array which has an element for each irrep of the
**               point group of the molecule (n.b. not just the ones
**               with a non-zero number of basis functions). each 
**               element contains the number of symmetry orbitals for
**               that irrep. Also, see chkpt_rd_orbspi().
**
** returns: none
**
** \ingroup (CHKPT)
*/

void chkpt_wt_sopi(int *sopi)
{
  int nirreps;
  char *keyword;
  keyword = chkpt_build_keyword("SO's per irrep");

  nirreps = chkpt_rd_nirreps();

  psio_write_entry(PSIF_CHKPT, keyword, (char *) sopi, nirreps*sizeof(int));

  free(keyword);
}
