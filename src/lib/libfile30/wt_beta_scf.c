#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_wt_beta_scf():  writes in the beta SCF eigenvector (or whatever is stored in
**     its place) to file30.
**
**   arguments: double **scf_vector    This rectangular matrix has dimentions nso
**     by nmo (see: rd_nmo()). For STO water, scf_vector 
**     should look something like the following:
**
**         *** *** *** *** 0.0 0.0 0.0
**         *** *** *** *** 0.0 0.0 0.0
**         *** *** *** *** 0.0 0.0 0.0
**         *** *** *** *** 0.0 0.0 0.0
**         0.0 0.0 0.0 0.0 *** 0.0 0.0
**         0.0 0.0 0.0 0.0 0.0 *** ***
**         0.0 0.0 0.0 0.0 0.0 *** ***
**
**     where the *** represent the non-zero values, and the 0.0 entries
**     represent (double)0.
**
**   returns: nothing.
*/

void file30_wt_beta_scf(double **scf_vector)
{
  int i,j,k;
  int mo_offset, so_offset;
  int nmo, mxcoef, nirreps, count;
  int *mopi, *sopi;
  PSI_FPTR junk, mo_coeff_ptr, scf_ptr;
  int *scf_ptrs;
  int tmp;
  double *wt_vector;

  /* STB(11/1/99) - Added to utilize the pointer array now in the SCF section
     of file30 */
  scf_ptrs = file30_rd_scf_ptrs();
  mo_coeff_ptr = (PSI_FPTR) (scf_ptrs[1] -1)*sizeof(int);
  
  /*scf_ptr = (PSI_FPTR) (info30_.mcalcs[0] + 60 - 1)*sizeof(int);
  wreadw(info30_.filenum, (char *) &tmp, sizeof(int), scf_ptr, &junk);
  mo_coeff_ptr = (PSI_FPTR) (tmp - 1)*sizeof(int);*/

  nirreps = file30_rd_nirreps();
  mxcoef = file30_rd_mxcoef();
  nmo = file30_rd_nmo();
  mopi = file30_rd_orbspi();
  sopi = file30_rd_sopi();
  wt_vector = init_array(mxcoef);

  count = 0;
  mo_offset = 0; so_offset = 0;
  for(i=0; i < nirreps; i++) 
   {
    for(j=mo_offset;j<(mo_offset+mopi[i]);j++)
     {
      for(k=so_offset;k<so_offset+sopi[i];k++,count++)
        {
         wt_vector[count] = scf_vector[k][j];
        }
     }
    so_offset += sopi[i];
    mo_offset += mopi[i];
   }

  wwritw(info30_.filenum, (char *) wt_vector, mxcoef*sizeof(double),
         mo_coeff_ptr, &junk);

  free(wt_vector); free(mopi); free(sopi); 
} 
  
  

