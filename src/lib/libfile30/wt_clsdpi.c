#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/*
** file30_wt_clsdpi():  Writes out number of closed shells per irrep.
**
**   arguments:
**     int *clsdpi -- an array containing the number of closed shells for
**    every irrep in the point group, not just those used in the HF procedure.
**
**   returns: nothing
**     
*/

void file30_wt_clsdpi(int *clsdpi)
{
  int i,j;
  int nsymhf, nirreps, nmo, mxcoef;
  int *index;
  char **irr_labs;
  char **hfsym_labs;
  PSI_FPTR junk;
  PSI_FPTR clsdpi_ptr, scf_ptr;
  int *scf_ptrs;
  int tmp;

  nsymhf = file30_rd_nsymhf();
  nirreps = file30_rd_nirreps();
  mxcoef = file30_rd_mxcoef();
  nmo = file30_rd_nmo();
  irr_labs = file30_rd_irr_labs();
  hfsym_labs = file30_rd_hfsym_labs();
  scf_ptrs = file30_rd_scf_ptrs();
  index = init_int_array(nsymhf);

  /*scf_ptr = (PSI_FPTR)(info30_.mcalcs[0]+60-1)*sizeof(int);
  wreadw(info30_.filenum, (char *) &tmp, sizeof(int), scf_ptr,
     &junk);
  clsdpi_ptr = (PSI_FPTR) (tmp - 1)*sizeof(int);
  clsdpi_ptr += (PSI_FPTR) (mxcoef+nmo)*sizeof(double) + (PSI_FPTR) 2*nsymhf*sizeof(int);*/
  
  /* STB(11/1/99) - Added to utilize the pointer array now in the SCF section
     of file30 */
  
  clsdpi_ptr = (PSI_FPTR) (scf_ptrs[6] - 1) *sizeof(int);

  for(i=0;i<nsymhf;i++)
    {index[i] = 0;
     for(j=i;j<nirreps;j++)
        if(!(strcmp(hfsym_labs[i],irr_labs[j]))) index[i] = j;
    }

  for(i=0;i<nsymhf;i++)
     wwritw(info30_.filenum, (char *) &clsdpi[index[i]], (int) sizeof(int),
	 clsdpi_ptr, &clsdpi_ptr);

  free(index);
  for(i=0; i < nirreps; i++) free(irr_labs[i]);
  free(irr_labs);
  for(i=0; i < nsymhf; i++) free(hfsym_labs[i]);
  free(hfsym_labs);
  free(scf_ptrs);

  return;
}
