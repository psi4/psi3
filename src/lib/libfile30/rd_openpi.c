/*!
  \file rd_openpi.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr/libciomr.h>

/*!
** file30_rd_openpi():  Reads in the number of open shells per irrep. 
**
**   takes no arguments.
**
**   returns:
**     int *openpi  an array which has an element for each irrep of the
**          point group of the molecule (n.b. not just the ones
**          with a non-zero number of basis functions).  Each element
**          contains the number of open shells for that irrep.
*/


int *file30_rd_openpi(void)
{
  int i, j;
  int nsymhf, nirreps, nmo, mxcoef;
  int *index, *openpi;
  char **irr_labs;
  char **hfsym_labs;
  PSI_FPTR junk;
  PSI_FPTR openpi_ptr, scf_ptr;
  int *scf_ptrs;
  int tmp;

  nsymhf = file30_rd_nsymhf();
  nirreps = file30_rd_nirreps();
  mxcoef = file30_rd_mxcoef();
  nmo = file30_rd_nmo();
  irr_labs = file30_rd_irr_labs();
  hfsym_labs = file30_rd_hfsym_labs();
  
/* STB(10/29/99) - Added to utilize the pointer array now in the SCF section
     of file30 */
  scf_ptrs = file30_rd_scf_ptrs();

  index = init_int_array(nsymhf);
  openpi = init_int_array(nirreps);

  /*scf_ptr = (PSI_FPTR) (info30_.mcalcs[0]+60-1)*sizeof(int);
  wreadw(info30_.filenum, (char *) &tmp, sizeof(int), scf_ptr,
     &junk);
  openpi_ptr = (PSI_FPTR) (tmp - 1)*sizeof(int);
  openpi_ptr += (mxcoef+nmo)*sizeof(double) +
  (PSI_FPTR) (nsymhf*4*sizeof(char) + 2*nsymhf*sizeof(int));*/

/* STB(10/29/99) - Added to utilize the pointer array now in the SCF section
     of file30 */
  openpi_ptr = (PSI_FPTR) (scf_ptrs[7] -1)*sizeof(int);

  for(i=0;i<nsymhf;i++)
    {index[i] = 0;
     for(j=i;j<nirreps;j++)
        if(!(strcmp(hfsym_labs[i],irr_labs[j]))) index[i] = j;
    }

  if(file30_rd_iopen()) {
     for(i=0;i<nsymhf;i++)
        wreadw(info30_.filenum, (char *) &openpi[index[i]], (int) sizeof(int),
               openpi_ptr, &openpi_ptr);
   }

  free(index);
  for(i=0; i < nirreps; i++) free(irr_labs[i]);
  free(irr_labs);
  for(i=0; i < nsymhf; i++) free(hfsym_labs[i]);
  free(hfsym_labs);
  free(scf_ptrs);

  return openpi;
}
