/*!
  \file cdsalc2cd.c
  \ingroup (CHKPT)
*/

#include <stdio.h>
#include <stdlib.h>
#include <psifiles.h>
#include "chkpt.h"
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_cdsalc2cd(): Read in (normalized) SALCs of cartesian displacements
**
** takes no arguments.
**
** returns: cdsalc2cd = A natom*3 by natom*3 blocked matrix of doubles. columnts correpond to symmetry-blocked SALCs
** 
** \ingroup (CHKPT)
*/

double **chkpt_rd_cdsalc2cd(void)
{
  const int num_cd = 3*chkpt_rd_natom();
  double **cdsalc2cd = block_matrix(num_cd,num_cd);
  psio_address ptr = PSIO_ZERO;
  char *keyword = chkpt_build_keyword("cartdisp->SALC matrix");

  psio_read(PSIF_CHKPT, keyword, (char *) cdsalc2cd[0], num_cd*num_cd*sizeof(double), ptr, &ptr);

  free(keyword);
  return cdsalc2cd;
}


/*!
** chkpt_wt_cdsalc2cd(): Writes out (normalized) SALCs of cartesian displacements
**
** \param cdsalc2cd = A natom*3 by natom*3 blocked matrix of doubles. columnts correpond to symmetry-blocked SALCs
**
** returns: none
**
** \ingroup (CHKPT)
*/

void chkpt_wt_cdsalc2cd(const double **cdsalc2cd)
{
  const int num_cd = 3*chkpt_rd_natom();
  psio_address ptr = PSIO_ZERO;
  char *keyword = chkpt_build_keyword("cartdisp->SALC matrix");

  psio_write(PSIF_CHKPT, keyword, (char *) cdsalc2cd[0], num_cd*num_cd*sizeof(double), ptr, &ptr);

  free(keyword);
}

