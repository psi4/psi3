/*!
  \file schwartz.c
  \ingroup (CHKPT)
*/

#include <stdio.h>
#include "chkpt.h"
#include <libciomr/libciomr.h>
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_schwartz(): Read in the table of maxima of Schwartz integrals 
**                      (ij|ij) for each shell doublet;
**
**  takes no arguments.
**
**  returns: NULL if no table is present in the file, double** (num_shells by 
**           num_shells) otherwise
**
**               ** MAY THE SCHWARTZ BE WITH YOU!!! **
**
** \ingroup (CHKPT)
*/
double **chkpt_rd_schwartz(void)
{
  double **schwartz;
  int num_shells;

  if (!psio_tocscan(PSIF_CHKPT, "::Schwartz table"))
    return NULL;
  else {
    num_shells = file30_rd_nshell();
    schwartz = block_matrix(num_shells, num_shells);

    psio_read_entry(PSIF_CHKPT, "::Schwartz table", (char *) schwartz[0],
		    num_shells*num_shells*sizeof(double));

    return schwartz;
  }
  
}


/*!
** chkpt_wt_schwartz(): Write out the table of maxima of Schwartz integrals 
**                      (ij|ij) for each shell doublet;
**
** \param schwartz = matrix (num_shells by num_shells)
**
**               ** MAY THE SCHWARTZ BE WITH YOU!!! **
**
** \ingroup (CHKPT)
*/
void chkpt_wt_schwartz(double **schwartz)
{
  int num_shells;

  num_shells = file30_rd_nshell();

  psio_write_entry(PSIF_CHKPT, "::Schwartz table", (char *) schwartz[0],
		  num_shells*num_shells*sizeof(double));
}

