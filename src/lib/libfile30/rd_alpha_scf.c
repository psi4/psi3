/*!
  \file rd_alpha_scf.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/*!
** file30_rd_alpha_scf():  Reads in the alpha SCF eigenvector (or whatever is stored in its
**   place).
**  
**   takes no arguments.
**  
**   returns: double **scf_vector    This rectangular matrix has dimentions nso
**     by nmo (see: rd_nmo()). For STO water, scf_vector would 
**     come out looking something like the following:
**        
**         *** *** *** *** 0.0 0.0 0.0
**         *** *** *** *** 0.0 0.0 0.0
**         *** *** *** *** 0.0 0.0 0.0
**         *** *** *** *** 0.0 0.0 0.0
**         0.0 0.0 0.0 0.0 *** 0.0 0.0
**         0.0 0.0 0.0 0.0 0.0 *** ***
**         0.0 0.0 0.0 0.0 0.0 *** ***
**
**    where the *** represent the non-zero values, and the 0.0 entries 
**    represent (double)0.
*/


double **file30_rd_alpha_scf(void)
{
  int i,j;
  int irrep, so, mo, so_offset, mo_offset;
  double **scf_vector;
  int nmo, mxcoef, nirreps, count, num_orbs;
  int sym, first, last, column;
  int *orbsym, *mopi, *sopi, *start;
  PSI_FPTR junk, mo_coeff_ptr;
  PSI_FPTR scf_ptr;
  int *scf_ptrs;
  int tmp;
  double *tmp_vector;

  /*scf_ptr = (PSI_FPTR) (info30_.mcalcs[0] + 60 - 1)*sizeof(int);*/
  /*wreadw(info30_.filenum, (char *) &tmp, sizeof(int), scf_ptr, &junk);*/
  scf_ptrs = file30_rd_scf_ptrs();
  mo_coeff_ptr = (PSI_FPTR) (scf_ptrs[0] -1)*sizeof(int);
  
/*mo_coeff_ptr = (PSI_FPTR) (tmp - 1)*sizeof(int);*/

  nirreps = file30_rd_nirreps();
  mxcoef = file30_rd_mxcoef();
  tmp_vector = init_array(mxcoef);
  wreadw(info30_.filenum, (char *) tmp_vector, mxcoef*sizeof(double),
	 mo_coeff_ptr, &junk);

  mopi = file30_rd_orbspi();
  sopi = file30_rd_sopi();
  scf_vector = block_matrix(file30_rd_nso(),file30_rd_nmo());

/*  count = 0;
  for(i=0; i < nirreps; i++) {
      start[i] = count;
      count += opi[i];
    }
  count = -1;
  for(i=0; i < nirreps; i++) {
      num_orbs = opi[i];
      for(j=0; j < num_orbs; j++) {
          count++;
          orbsym[count] = i;
        }
    }
  count = -1;
  for(i=0; i < nmo; i++) {
      column = i;
      sym = orbsym[column];
      first = start[sym];
      last = start[sym] + opi[sym];
      for(j=first; j < last; j++) {
          count += 1;
          scf_vector[j][column] = tmp_vector[count];
        }
    }*/

  count = 0;
  so_offset = 0; mo_offset = 0;
  for(irrep=0;irrep < nirreps; irrep++)
    if (sopi[irrep] > 0) {
      for(mo=0; mo<mopi[irrep]; mo++)
	for(so=0; so<sopi[irrep]; so++) {
	  scf_vector[so+so_offset][mo+mo_offset] = tmp_vector[count];
	  count++;
	}
      so_offset += sopi[irrep];
      mo_offset += mopi[irrep];
    }
      

  free(sopi);  free(mopi); free(tmp_vector);
  free(scf_ptrs);

  return scf_vector;
}
