/*!
  \file wt_openpi.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/*!
** file30_wt_openpi():  Writes out number of open shells per irrep.
**
**   arguments:
** \param int *openpi -- an array containing the number of open shells for
**    every irrep in the point group, not just those used in the HF procedure.
**
**   returns: nothing
**     
*/

void file30_wt_openpi(int *openpi)
{
  int i,j;
  int nsymhf, nirreps, nmo, mxcoef;
  int *index;
  char **irr_labs;
  char **hfsym_labs;
  int *scf_ptrs;
  PSI_FPTR junk;
  PSI_FPTR openpi_ptr, scf_ptr;
  int tmp;

  nsymhf = file30_rd_nsymhf();
  nirreps = file30_rd_nirreps();
  mxcoef = file30_rd_mxcoef();
  nmo = file30_rd_nmo();
  irr_labs = file30_rd_irr_labs();
  scf_ptrs = file30_rd_scf_ptrs();
  hfsym_labs = file30_rd_hfsym_labs();

  index = init_int_array(nsymhf);

  /*scf_ptr = (PSI_FPTR)(info30_.mcalcs[0]+60-1)*sizeof(int);
  wreadw(info30_.filenum, (char *) &tmp, sizeof(int), scf_ptr,
     &junk);
  openpi_ptr = (PSI_FPTR) (tmp - 1)*sizeof(int);
  openpi_ptr += (PSI_FPTR) (mxcoef+nmo)*sizeof(double) + (PSI_FPTR) 3*nsymhf*sizeof(int);*/
  /* STB(11/1/99) - Added to utilize the pointer array now in the SCF section
     of file30 */

  openpi_ptr = (PSI_FPTR) (scf_ptrs[7] - 1)*sizeof(int);

  for(i=0;i<nsymhf;i++)
    {index[i] = 0;
     for(j=i;j<nirreps;j++)
        if(!(strcmp(hfsym_labs[i],irr_labs[j]))) index[i] = j;
    }

  for(i=0;i<nsymhf;i++)
     wwritw(info30_.filenum, (char *) &openpi[index[i]], (int) sizeof(int),
	 openpi_ptr, &openpi_ptr);

  free(index);
  for(i=0; i < nirreps; i++) free(irr_labs[i]);
  free(irr_labs);
  for(i=0; i < nsymhf; i++) free(hfsym_labs[i]);
  free(hfsym_labs);
  free(scf_ptrs);

  return;
}
