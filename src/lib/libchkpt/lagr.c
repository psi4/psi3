/*!
  \file lagr.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "chkpt.h"
#include <libciomr/libciomr.h>
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_lagr():  Reads in the MO lagrangian matrix for RHF/ROHF.
**
**   takes no arguments.
**
**   returns: 
**	double **lagr	a matrix nmo by nmo.
*/
double **chkpt_rd_lagr(void)
{
  int nmo;
  double **lagr;

  nmo = chkpt_rd_nmo();

  lagr = block_matrix(nmo,nmo);
  psio_read_entry(PSIF_CHKPT, "::MO Lagrangian", (char *) lagr[0], 
		  nmo*nmo*sizeof(double));

  return lagr;
}

/*!
** chkpt_wt_lagr():  Writes the MO lagrangian matrix for RHF/ROHF.
**
**  argument: 
**    \param double **lagr	a matrix nmo by nmo.
**
** returns: none
*/
void chkpt_wt_lagr(double **lagr)
{
  int nmo;

  nmo = chkpt_rd_nmo();

  psio_write_entry(PSIF_CHKPT, "::MO Lagrangian", (char *) lagr[0], 
		   nmo*nmo*sizeof(double));
}


/*!
** chkpt_rd_alpha_lagr():  Reads in the alpha MO lagrangian matrix for UHF.
**
**   takes no arguments.
**
**   returns: 
**	double **lagr	a matrix nmo by nmo.
*/
double **chkpt_rd_alpha_lagr(void)
{
  int nmo;
  double **lagr;

  nmo = chkpt_rd_nmo();

  lagr = block_matrix(nmo,nmo);
  psio_read_entry(PSIF_CHKPT, "::Alpha MO Lagrangian", (char *) lagr[0], 
		  nmo*nmo*sizeof(double));

  return lagr;
}

/*!
** chkpt_wt_alpha_lagr():  Writes the alpha MO lagrangian matrix for UHF.
**
**  argument: 
**    \param double **lagr	a matrix nmo by nmo.
**
** returns: none
*/
void chkpt_wt_alpha_lagr(double **lagr)
{
  int nmo;

  nmo = chkpt_rd_nmo();

  psio_write_entry(PSIF_CHKPT, "::Alpha MO Lagrangian", (char *) lagr[0], 
		   nmo*nmo*sizeof(double));
}

/*!
** chkpt_rd_beta_lagr():  Reads in the beta MO lagrangian matrix for UHF.
**
**   takes no arguments.
**
**   returns: 
**	double **lagr	a matrix nmo by nmo.
*/
double **chkpt_rd_beta_lagr(void)
{
  int nmo;
  double **lagr;

  nmo = chkpt_rd_nmo();

  lagr = block_matrix(nmo,nmo);
  psio_read_entry(PSIF_CHKPT, "::Beta MO Lagrangian", (char *) lagr[0], 
		  nmo*nmo*sizeof(double));

  return lagr;
}

/*!
** chkpt_wt_beta_lagr():  Writes the beta MO lagrangian matrix for UHF.
**
**  argument: 
**    \param double **lagr	a matrix nmo by nmo.
**
** returns: none
*/
void chkpt_wt_beta_lagr(double **lagr)
{
  int nmo;

  nmo = chkpt_rd_nmo();

  psio_write_entry(PSIF_CHKPT, "::Beta MO Lagrangian", (char *) lagr[0], 
		   nmo*nmo*sizeof(double));
}
