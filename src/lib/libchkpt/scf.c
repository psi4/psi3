/*!
  \file scf.c
*/

#include <stdio.h>
#include "chkpt.h"
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <psifiles.h>

/*!
** chkpt_rd_scf():  Reads in the full SCF eigenvector matrix for RHF/ROHF.
**  
**   takes no arguments.
**  
**   returns: double **scf    This rectangular matrix has dimensions nso
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
  double **scf;
  int nmo, nso;

  nmo = chkpt_rd_nmo();
  nso = chkpt_rd_nso();

  scf = block_matrix(nso,nmo);
  psio_read_entry(PSIF_CHKPT, "::MO coefficients", (char *) scf[0], nso*nmo*sizeof(double));

  return scf;
}

/*!
** chkpt_rd_alpha_scf():  Reads in the full alpha SCF eigenvector matrix for UHF.
**  
**   takes no arguments.
**  
**   returns: double **scf    This rectangular matrix has dimensions nso
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
  double **scf;
  int nmo, nso;

  nmo = chkpt_rd_nmo();
  nso = chkpt_rd_nso();

  scf = block_matrix(nso,nmo);
  psio_read_entry(PSIF_CHKPT, "::Alpha MO coefficients", (char *) scf[0], nso*nmo*sizeof(double));

  return scf;
}

/*!
** chkpt_rd_beta_scf():  Reads in the full beta SCF eigenvector matrix for UHF.
**  
**   takes no arguments.
**  
**   returns: double **scf    This rectangular matrix has dimensions nso
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
  double **scf;
  int nmo, nso;

  nmo = chkpt_rd_nmo();
  nso = chkpt_rd_nso();

  scf = block_matrix(nso,nmo);
  psio_read_entry(PSIF_CHKPT, "::Beta MO coefficients", (char *) scf[0], nso*nmo*sizeof(double));

  return scf;
}

/*!
** chkpt_wt_scf():  Writes the full SCF eigenvector matrix for RHF/ROHF.
**  
** arguments: 
**  \param double **scf    This rectangular matrix has dimensions nso
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
**
** returns: none
**
** NOTE: The input scf matrix must occupy a contiguous block of nmo x
** nso memory.  Use block_matrix() from libciomr to allocate space for
** the matrix.
*/
void chkpt_wt_scf(double **scf)
{
  int nmo, nso;

  nmo = chkpt_rd_nmo();
  nso = chkpt_rd_nso();

  psio_write_entry(PSIF_CHKPT, "::MO coefficients", (char *) scf[0], nso*nmo*sizeof(double));
}

/*!
** chkpt_wt_alpha_scf():  Writes the full alpha SCF eigenvector matrix for UHF.
**  
** arguments: 
**  \param double **scf    This rectangular matrix has dimensions nso
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
**
** returns: none
**
** NOTE: The input scf matrix must occupy a contiguous block of nmo x
** nso memory.  Use block_matrix() from libciomr to allocate space for
** the matrix.
*/
void chkpt_wt_alpha_scf(double **scf)
{
  int nmo, nso;

  nmo = chkpt_rd_nmo();
  nso = chkpt_rd_nso();

  psio_write_entry(PSIF_CHKPT, "::Alpha MO coefficients", (char *) scf[0], nso*nmo*sizeof(double));
}

/*!
** chkpt_wt_beta_scf():  Writes the full beta SCF eigenvector matrix for UHF.
**  
** arguments: 
**  \param double **scf    This rectangular matrix has dimensions nso
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
**
** returns: none
**
** NOTE: The input scf matrix must occupy a contiguous block of nmo x
** nso memory.  Use block_matrix() from libciomr to allocate space for
** the matrix.
*/
void chkpt_wt_beta_scf(double **scf)
{
  int nmo, nso;

  nmo = chkpt_rd_nmo();
  nso = chkpt_rd_nso();

  psio_write_entry(PSIF_CHKPT, "::Beta MO coefficients", (char *) scf[0], nso*nmo*sizeof(double));
}


/*!
** chkpt_rd_scf_irrep(): Reads a single irrep of the SCF eigenvectors for RHF/ROHF.
**
** arguments:
**  \param int irrep  The desired irreducible representation.
**
** returns: double **scf   A rectangualr sopi[irrep] by orbspi[irrep] matrix.
*/

double **chkpt_rd_scf_irrep(int irrep)
{
  int i, j, row, col;
  int nirreps, nso, nmo;
  int *sopi, *mopi;
  double **scf, **scf_full;

  nirreps = chkpt_rd_nirreps();
  sopi = chkpt_rd_sopi();
  mopi = chkpt_rd_orbspi();
  nso = chkpt_rd_nso();
  nmo = chkpt_rd_nmo();

  scf = block_matrix(sopi[irrep],mopi[irrep]);
  scf_full = chkpt_rd_scf();

  /* compute row and column offsets */
  for(i=0,row=0,col=0; i < irrep; i++) {
    row += sopi[i];
    col += mopi[i];
  }

  for(i=0; i < sopi[irrep]; i++)
    for(j=0; j < mopi[irrep]; j++)
      scf[i][j] = scf_full[i+row][j+col];

  free_block(scf_full);
  free(sopi);
  free(mopi);

  return scf;
}

/*!
** chkpt_rd_scf_alpha_irrep(): Reads a single irrep of the alpha SCF eigenvectors for UHF.
**
** arguments:
**  \param int irrep  The desired irreducible representation.
**
** returns: double **scf   A rectangualr sopi[irrep] by orbspi[irrep] matrix.
*/

double **chkpt_rd_alpha_scf_irrep(int irrep)
{
  int i, j, row, col;
  int nirreps, nso, nmo;
  int *sopi, *mopi;
  double **scf, **scf_full;

  nirreps = chkpt_rd_nirreps();
  sopi = chkpt_rd_sopi();
  mopi = chkpt_rd_orbspi();
  nso = chkpt_rd_nso();
  nmo = chkpt_rd_nmo();

  scf = block_matrix(sopi[irrep],mopi[irrep]);
  scf_full = chkpt_rd_alpha_scf();

  /* compute row and column offsets */
  for(i=0,row=0,col=0; i < irrep; i++) {
    row += sopi[i];
    col += mopi[i];
  }

  for(i=0; i < sopi[irrep]; i++)
    for(j=0; j < mopi[irrep]; j++)
      scf[i][j] = scf_full[i+row][j+col];

  free_block(scf_full);
  free(sopi);
  free(mopi);

  return scf;
}

/*!
** chkpt_rd_scf_beta_irrep(): Reads a single irrep of the beta SCF eigenvectors for UHF.
**
** arguments:
**  \param int irrep  The desired irreducible representation.
**
** returns: double **scf   A rectangualr sopi[irrep] by orbspi[irrep] matrix.
*/

double **chkpt_rd_beta_scf_irrep(int irrep)
{
  int i, j, row, col;
  int nirreps, nso, nmo;
  int *sopi, *mopi;
  double **scf, **scf_full;

  nirreps = chkpt_rd_nirreps();
  sopi = chkpt_rd_sopi();
  mopi = chkpt_rd_orbspi();
  nso = chkpt_rd_nso();
  nmo = chkpt_rd_nmo();

  scf = block_matrix(sopi[irrep],mopi[irrep]);
  scf_full = chkpt_rd_beta_scf();

  /* compute row and column offsets */
  for(i=0,row=0,col=0; i < irrep; i++) {
    row += sopi[i];
    col += mopi[i];
  }

  for(i=0; i < sopi[irrep]; i++)
    for(j=0; j < mopi[irrep]; j++)
      scf[i][j] = scf_full[i+row][j+col];

  free_block(scf_full);
  free(sopi);
  free(mopi);

  return scf;
}

/*!
** chkpt_wt_scf_irrep(): Writes a single irrep of the SCF eigenvectors for RHF/ROHF.
**
** arguments:
**  \param int irrep  The desired irreducible representation.
**
** returns: double **scf   A rectangualr sopi[irrep] by orbspi[irrep] matrix.
*/

void chkpt_wt_scf_irrep(double **scf, int irrep)
{
  int i, j, row, col;
  int nirreps, nso, nmo;
  int *sopi, *mopi;
  double **scf_full;

  nirreps = chkpt_rd_nirreps();
  sopi = chkpt_rd_sopi();
  mopi = chkpt_rd_orbspi();
  nso = chkpt_rd_nso();
  nmo = chkpt_rd_nmo();

  scf_full = chkpt_rd_scf();

  /* compute row and column offsets */
  for(i=0,row=0,col=0; i < irrep; i++) {
    row += sopi[i];
    col += mopi[i];
  }

  for(i=0; i < sopi[irrep]; i++)
    for(j=0; j < mopi[irrep]; j++)
      scf_full[i+row][j+col] = scf[i][j];

  chkpt_wt_scf(scf_full);
  free_block(scf_full);
  free(sopi);
  free(mopi);
}

/*!
** chkpt_wt_alpha_scf_irrep(): Writes a single irrep of the alpha SCF 
**                             eigenvectors for RHF/ROHF.
**
** arguments:
**  \param int irrep  The desired irreducible representation.
**
** returns: double **scf   A rectangualr sopi[irrep] by orbspi[irrep] matrix.
*/
void chkpt_wt_alpha_scf_irrep(double **scf, int irrep)
{
  int i, j, row, col;
  int nirreps, nso, nmo;
  int *sopi, *mopi;
  double **scf_full;

  nirreps = chkpt_rd_nirreps();
  sopi = chkpt_rd_sopi();
  mopi = chkpt_rd_orbspi();
  nso = chkpt_rd_nso();
  nmo = chkpt_rd_nmo();

  scf_full = chkpt_rd_alpha_scf();

  /* compute row and column offsets */
  for(i=0,row=0,col=0; i < irrep; i++) {
    row += sopi[i];
    col += mopi[i];
  }

  for(i=0; i < sopi[irrep]; i++)
    for(j=0; j < mopi[irrep]; j++)
      scf_full[i+row][j+col] = scf[i][j];

  chkpt_wt_alpha_scf(scf_full);
  free_block(scf_full);
  free(sopi);
  free(mopi);
}

/*!
** chkpt_wt_beta_scf_irrep(): Writes a single irrep of the beta SCF 
**                             eigenvectors for RHF/ROHF.
**
** arguments:
**  \param int irrep  The desired irreducible representation.
**
** returns: double **scf   A rectangualr sopi[irrep] by orbspi[irrep] matrix.
*/
void chkpt_wt_beta_scf_irrep(double **scf, int irrep)
{
  int i, j, row, col;
  int nirreps, nso, nmo;
  int *sopi, *mopi;
  double **scf_full;

  nirreps = chkpt_rd_nirreps();
  sopi = chkpt_rd_sopi();
  mopi = chkpt_rd_orbspi();
  nso = chkpt_rd_nso();
  nmo = chkpt_rd_nmo();

  scf_full = chkpt_rd_beta_scf();

  /* compute row and column offsets */
  for(i=0,row=0,col=0; i < irrep; i++) {
    row += sopi[i];
    col += mopi[i];
  }

  for(i=0; i < sopi[irrep]; i++)
    for(j=0; j < mopi[irrep]; j++)
      scf_full[i+row][j+col] = scf[i][j];

  chkpt_wt_beta_scf(scf_full);
  free_block(scf_full);
  free(sopi);
  free(mopi);
}


