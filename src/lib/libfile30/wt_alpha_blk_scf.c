/*!
  \file wt_alpha_blk_scf.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr/libciomr.h>

/*!
** file30_wt_alpha_blk_scf():  Writes in the alpha SCF eigenvector (or whateve 
**     is to be stored in its place).
**
**   arguments: 
**
** \param int irrep   the number of the irrep to which the symmetry block 
**        belongs (this includes irreps with orbspi[irrep] == 0)
**        n.b. this routine assumes that the first irrep will have
**        irrep == 0.
**
** \param double **scf_vector    This should be a single symmetry
**        block of the SCF eigenvector.  Its dimension should be 
**        sopi[irrep]*orbspi[irrep];
**
**   returns: nothing.
*/

void file30_wt_alpha_blk_scf(double **scf_vector, int irrep)
{
  int i,j,k;
  int nmo, mxcoef, nirreps, count, offset;
  int *mopi, *sopi;
  PSI_FPTR junk, mo_coeff_ptr, scf_ptr;
  int *scf_ptrs;
  int tmp;
  double *wt_vector;

  /*scf_ptr = (PSI_FPTR) (info30_.mcalcs[0] + 60 - 1)*sizeof(int);
  wreadw(info30_.filenum, (char *) &tmp, sizeof(int), scf_ptr, &junk);
  mo_coeff_ptr = (PSI_FPTR) (tmp - 1)*sizeof(int);*/

/* STB(11/1/99) - Added to utilize the pointer array now in the SCF section
     of file30 */
  scf_ptrs = file30_rd_scf_ptrs();
  mo_coeff_ptr = (PSI_FPTR) (scf_ptrs[0]-1)*sizeof(int);

  mopi = file30_rd_orbspi();
  sopi = file30_rd_sopi();
  
  if(sopi[irrep]) 
   {
    offset = 0;
    for(i=0; i < irrep; i++) {
        offset += sopi[i]*mopi[i];
      } 
    mo_coeff_ptr += (PSI_FPTR) offset * sizeof(double);

    count = 0;
    wt_vector = init_array(sopi[irrep]*mopi[irrep]);
    for(j=0;j<mopi[irrep];j++)
       {
        for(k=0;k<sopi[irrep];k++,count++)
          {
           wt_vector[count] = scf_vector[k][j];
          }
       }
  
    wwritw(info30_.filenum, (char *) wt_vector,
          sopi[irrep]*mopi[irrep]*sizeof(double), mo_coeff_ptr, &junk);

    free(wt_vector);
   }

  free(mopi);
  free(sopi);
  free(scf_ptrs);

  return;
} 
  
  

