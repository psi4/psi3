/*!
  \file evals.c
*/

#include <stdio.h>
#include <libciomr/libciomr.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_evals():  Reads in the SCF orbital energies for RHF/ROHF.
**
**  takes no arguments.
**
**  returns: double *evals   an array of _all_ of the SCF eigenvalues,
**      ordered by irrep, and by increasing energy within each irrep.  
**      (i.e. for sto water, the four a1 eigenvalues all come first, and 
**      those four are ordered from lowest energy to highest energy,
**      followed by the single b1 eigenvalue, etc.)
*/

double *chkpt_rd_evals(void)
{
  double *energies;
  
  energies = init_array(chkpt_rd_nmo());
  psio_read_entry(PSIF_CHKPT, "::MO energies", (char *) energies, 
		  chkpt_rd_nmo()*sizeof(double));

  return energies;
}


/*!
** chkpt_rd_alpha_evals():  Reads in the SCF alpha orbital energies for UHF.
**
**  takes no arguments.
**
**  returns: double *evals   an array of _all_ of the alpha SCF eigenvalues,
**      ordered by irrep, and by increasing energy within each irrep.  
**      (i.e. for sto water, the four a1 eigenvalues all come first, and 
**      those four are ordered from lowest energy to highest energy,
**      followed by the single b1 eigenvalue, etc.)
*/

double *chkpt_rd_alpha_evals(void)
{
  double *energies;
  
  energies = init_array(chkpt_rd_nmo());
  psio_read_entry(PSIF_CHKPT, "::Alpha MO energies", (char *) energies, 
		  chkpt_rd_nmo()*sizeof(double));

  return energies;
}

/*!
** chkpt_rd_beta_evals():  Reads in the SCF beta orbital energies for UHF.
**
**  takes no arguments.
**
**  returns: double *evals   an array of _all_ of the beta SCF eigenvalues,
**      ordered by irrep, and by increasing energy within each irrep.  
**      (i.e. for sto water, the four a1 eigenvalues all come first, and 
**      those four are ordered from lowest energy to highest energy,
**      followed by the single b1 eigenvalue, etc.)
*/

double *chkpt_rd_beta_evals(void)
{
  double *energies;
  
  energies = init_array(chkpt_rd_nmo());
  psio_read_entry(PSIF_CHKPT, "::Beta MO energies", (char *) energies, 
		  chkpt_rd_nmo()*sizeof(double));

  return energies;
}

/*!
** chkpt_wt_evals():  Writes the SCF orbital energies for UHF.
**
** arguments: 
**  \param double *evals  an array of _all_ of the SCF eigenvalues,
**      ordered by irrep, and by increasing energy within each irrep.  
**      (i.e. for sto water, the four a1 eigenvalues all come first, and 
**      those four are ordered from lowest energy to highest energy,
**      followed by the single b1 eigenvalue, etc.)
**
** returns: none
*/

void chkpt_wt_evals(double *energies)
{
  psio_write_entry(PSIF_CHKPT, "::MO energies", (char *) energies, 
		   chkpt_rd_nmo()*sizeof(double));
}

/*!
** chkpt_wt_alpha_evals():  Writes the SCF alpha orbital energies for UHF.
**
** arguments: 
**  \param double *evals  an array of _all_ of the alpha SCF eigenvalues,
**      ordered by irrep, and by increasing energy within each irrep.  
**      (i.e. for sto water, the four a1 eigenvalues all come first, and 
**      those four are ordered from lowest energy to highest energy,
**      followed by the single b1 eigenvalue, etc.)
**
** returns: none
*/

void chkpt_wt_alpha_evals(double *energies)
{
  psio_write_entry(PSIF_CHKPT, "::Alpha MO energies", (char *) energies, 
		   chkpt_rd_nmo()*sizeof(double));
}

/*!
** chkpt_wt_beta_evals():  Writes the SCF beta orbital energies for UHF.
**
** arguments: 
**  \param double *evals  an array of _all_ of the beta SCF eigenvalues,
**      ordered by irrep, and by increasing energy within each irrep.  
**      (i.e. for sto water, the four a1 eigenvalues all come first, and 
**      those four are ordered from lowest energy to highest energy,
**      followed by the single b1 eigenvalue, etc.)
**
** returns: none
*/

void chkpt_wt_beta_evals(double *energies)
{
  psio_write_entry(PSIF_CHKPT, "::Beta MO energies", (char *) energies, 
		   chkpt_rd_nmo()*sizeof(double));
}

