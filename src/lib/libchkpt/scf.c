/*!
  \file scf.c
*/

#include <stdio.h>
#include "chkpt.h"
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <psifiles.h>

/*!
** chkpt_rd_scf_packed(): Reads the symmetry-packed SCF vector for RHF/ROHF.
**
** takes no argument.
**
** returns: double *scf_vector  
**
*/
double *chkpt_rd_scf_packed(void)
{
  int mxcoef;
  double *scf_vector;

  mxcoef = chkpt_rd_maxcoef();
  scf_vector = init_array(maxcoef);

  psio_read_entry(PSIF_CHKPT, "::MO's packed", (char *) scf_vector, sizeof(double)*mxcoef);

  return scf_vector;
}

/*!
** chkpt_wt_scf_packed(): Writes the symmetry-packed SCF vector for RHF/ROHF.
**
** arguments: 
** \param double *scf_vector  
**
** returns: none
*/
void chkpt_wt_scf_packed(double *scf_vector)
{
  int mxcoef;

  mxcoef = chkpt_rd_maxcoef();

  psio_write_entry(PSIF_CHKPT, "::MO's packed", (char *) scf_vector, sizeof(double)*mxcoef);
}

/*!
** chkpt_rd_alpha_scf_packed(): Reads the symmetry-packed alpha SCF vector for UHF.
**
** takes no argument.
**
** returns: double *scf_vector  
**
*/
double *chkpt_rd_alpha_scf_packed(void)
{
  int mxcoef;
  double *scf_vector;

  mxcoef = chkpt_rd_maxcoef();
  scf_vector = init_array(maxcoef);

  psio_read_entry(PSIF_CHKPT, "::MO's alpha packed", (char *) scf_vector, 
		  sizeof(double)*mxcoef);

  return scf_vector;
}

/*!
** chkpt_wt_alpha_scf_packed(): Writes the symmetry-packed alpha SCF vector for UHF.
**
** arguments: 
** \param double *scf_vector  
**
** returns: none
*/
void chkpt_wt_alpha_scf_packed(double *scf_vector)
{
  int mxcoef;

  mxcoef = chkpt_rd_maxcoef();

  psio_write_entry(PSIF_CHKPT, "::MO's alpha packed", (char *) scf_vector, 
		   sizeof(double)*mxcoef);
}


/*!
** chkpt_rd_beta_scf_packed(): Reads the symmetry-packed beta SCF vector for UHF.
**
** takes no argument.
**
** returns: double *scf_vector  
**
*/
double *chkpt_rd_beta_scf_packed(void)
{
  int mxcoef;
  double *scf_vector;

  mxcoef = chkpt_rd_maxcoef();
  scf_vector = init_array(maxcoef);

  psio_read_entry(PSIF_CHKPT, "::MO's beta packed", (char *) scf_vector, 
		  sizeof(double)*mxcoef);

  return scf_vector;
}

/*!
** chkpt_wt_beta_scf_packed(): Writes the symmetry-packed beta SCF vector for UHF.
**
** arguments: 
** \param double *scf_vector  
**
** returns: none
*/
void chkpt_wt_beta_scf_packed(double *scf_vector)
{
  int mxcoef;

  mxcoef = chkpt_rd_maxcoef();

  psio_write_entry(PSIF_CHKPT, "::MO's beta packed", (char *) scf_vector, 
		   sizeof(double)*mxcoef);
}

/*!
** chkpt_rd_scf():  Reads in the full SCF eigenvector matrix for RHF/ROHF.
**  
**   takes no arguments.
**  
**   returns: double **scf_vector    This rectangular matrix has dimensions nso
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

double **chkpt_rd_scf(void)
{
  int irrep, so, mo, so_offset, mo_offset;
  double **scf_vector;
  int nmo, nirreps, count, num_orbs;
  int *mopi, *sopi;

  double *tmp_vector;

  nirreps = chkpt_rd_nirreps();
  tmp_vector = chkpt_rd_scf_packed();

  mopi = chkpt_rd_orbspi();
  sopi = chkpt_rd_sopi();
  scf_vector = block_matrix(file30_rd_nso(),file30_rd_nmo());

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

  return scf_vector;
}

/*!
** chkpt_rd_alpha_scf():  Reads in the alpha SCF eigenvectors for UHF.
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

double **chkpt_rd_alpha_scf(void)
{
  int irrep, so, mo, so_offset, mo_offset;
  double **scf_vector;
  int nmo, mxcoef, nirreps, count, num_orbs;
  int *mopi, *sopi;

  double *tmp_vector;

  nirreps = file30_rd_nirreps();
  mxcoef = file30_rd_mxcoef();
  tmp_vector = init_array(mxcoef);
  psio_read_entry(PSIF_CHKPT, "::Alpha MO coefficients", (char *) tmp_vector, 
		  sizeof(double)*mxcoef);

  mopi = chkpt_rd_orbspi();
  sopi = chkpt_rd_sopi();
  scf_vector = block_matrix(file30_rd_nso(),file30_rd_nmo());

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

  return scf_vector;
}

/*!
** chkpt_rd_beta_scf():  Reads in the beta SCF eigenvectors UHF.
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

double **chkpt_rd_beta_scf(void)
{
  int irrep, so, mo, so_offset, mo_offset;
  double **scf_vector;
  int nmo, mxcoef, nirreps, count, num_orbs;
  int *mopi, *sopi;

  double *tmp_vector;

  nirreps = file30_rd_nirreps();
  mxcoef = file30_rd_mxcoef();
  tmp_vector = init_array(mxcoef);
  psio_read_entry(PSIF_CHKPT, "::Beta MO coefficients", (char *) tmp_vector, 
		  sizeof(double)*mxcoef);

  mopi = chkpt_rd_orbspi();
  sopi = chkpt_rd_sopi();
  scf_vector = block_matrix(file30_rd_nso(),file30_rd_nmo());

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

  return scf_vector;
}

/*!
** chkpt_rd_block_scf():  Reads RHF SCF eigenvectors as a full matrix.
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

double **chkpt_rd_block_scf(void)
{
  int nmo, nso;
  double **scf_vector;
  psio_address ptr;

  nmo = chkpt_rd_nmo();
  nso = chkpt_rd_nso();

  ptr = PSIO_ZERO;
  for(i=0; i < nso; i++)
    psio_read(PSIF_CHKPT, "::MO coefficients (full block)", (char *) scf_vector[i],
	      sizeof(double)*nmo, ptr, &ptr);

  return scf_vector;
}

void chkpt_wt_block_scf(double **scf_vector)
{
  int i, nso, nmo;
  psio_address ptr;

  nmo = chkpt_rd_nmo();
  nso = chkpt_rd_nso();

  ptr = PSIO_ZERO;
  for(i=0; i < nso; i++)
    psio_write(PSIF_CHKPT, "::MO coefficients (full block)", (char *) scf_vector[i], 
	       sizeof(double)*nmo, ptr, &ptr);
}

