/*!
  \file rd_clsdpi.c
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
  int i, j;
  int nsymhf, nirreps;
  int *index, *clsdpi;
  char **irr_labs, **hfsym_labs;
  psio_address ptr;

  nsymhf = chkpt_rd_nsymhf();
  nirreps = chkpt_rd_nirreps();

  irr_labs = chkpt_rd_irr_labs();
  hfsym_labs = chkpt_rd_hfsym_labs();

  index = init_int_array(nsymhf);
  clsdpi = init_int_array(nirreps);

  for(i=0;i<nsymhf;i++)
    {index[i] = 0;
     for(j=i;j<nirreps;j++)
        if(!(strcmp(hfsym_labs[i],irr_labs[j]))) index[i] = j;
    }

  ptr = PSIO_ZERO;
  for(i=0;i<nsymhf;i++)
    psio_read(PSIF_CHKPT, "::Closed shells per HF irrep", (char *) clsdpi[index[i]], 
	      sizeof(int), ptr, &ptr);

  free(index);
  for(i=0; i < nirreps; i++) free(irr_labs[i]);
  free(irr_labs);
  for(i=0; i < nsymhf; i++) free(hfsym_labs[i]);
  free(hfsym_labs);
  
  return clsdpi;
}
