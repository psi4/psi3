#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_rd_beta_blk_scf():  Reads in a symmetry block of the beta SCF eigenvector 
**   (or whatever is stored in its place).
**  
**   arguments: int irrep   designates the particular irrep to which the block
**     belongs.
**  
**   returns: double **scf_vector    This rectangular matrix has orbspi[irrep] 
**     rows.
*/


double **file30_rd_beta_blk_scf(int irrep)
{
  int i,j;
  double **scf_vector;
  int nmo, mxcoef, nirreps, count, num_orbs;
  int sym, first, last, column, offset;
  int *mopi,*sopi;
  PSI_FPTR junk, mo_coeff_ptr, scf_ptr;
  int *scf_ptrs;
  int tmp;
  double *tmp_vector;

  /* STB(10/29/99) - Added to utilize the pointer array now in the SCF section
     of file30 */
  scf_ptrs = file30_rd_scf_ptrs();
  mo_coeff_ptr = (PSI_FPTR) (scf_ptrs[1]-1)*sizeof(int);
  /*scf_ptr = (PSI_FPTR) (info30_.mcalcs[0] + 60 - 1)*sizeof(int);
  wreadw(info30_.filenum, (char *) &tmp, sizeof(int), scf_ptr, &junk);
  mo_coeff_ptr = (PSI_FPTR) (tmp - 1)*sizeof(int);*/

  nirreps = file30_rd_nirreps();
  mxcoef = file30_rd_mxcoef();
  mopi = file30_rd_orbspi();
  sopi = file30_rd_sopi();

  scf_vector = NULL;

  if(sopi[irrep])
   {
    tmp_vector = init_array(sopi[irrep]*mopi[irrep]);

    offset = 0;
    for(i=0; i < irrep; i++) {
        offset += sopi[i]*mopi[i];
      }
    mo_coeff_ptr += (PSI_FPTR) offset*sizeof(double);
  
    wreadw(info30_.filenum, (char *) tmp_vector, 
           (sopi[irrep]*mopi[irrep])*sizeof(double), mo_coeff_ptr, &junk);
  
    scf_vector = init_matrix(sopi[irrep],mopi[irrep]);
    
    count = 0;
    for(i=0; i < mopi[irrep] ; i++)
       for(j=0; j < sopi[irrep] ; j++, count++) {
          scf_vector[j][i] = tmp_vector[count];
         }
    free(tmp_vector);
    }

  free(mopi);
  free(sopi);

  return scf_vector;
}


