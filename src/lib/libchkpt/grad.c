/*!
  \file grad.c
  \ingroup (CHKPT)
*/

#include <stdio.h>
#include <stdlib.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>

/*!
** chkpt_rd_grad():  Reads the energy gradient WRT nuclear coordinates
**
**   takes no arguments.
**
**   returns: grad = a vector of doubles natom*3 elements long, e.g.
**     grad[0] = gradient wrt x coordinate of atom 0
**     grad[1] = gradient wrt y coordinate of atom 0
**     grad[8] = gradient wrt z coordinate of atom 2
** \ingroup (CHKPT)
*/

double *chkpt_rd_grad(void)
{
  int natom;
  double *grad;
  char *keyword;
  keyword = chkpt_build_keyword("Energy Gradient");

  natom = chkpt_rd_natom();
  grad = init_array(natom*3);

  psio_read_entry(PSIF_CHKPT, keyword, (char *) grad, natom*3*sizeof(double));

  free(keyword);
  return grad;
}


/*!
** chkpt_wt_grad():  Writes the energy gradient WRT nuclear coordinates
**
**   arguments:
**   \param grad = a vector of doubles natom*3 elements long, e.g.
**     grad[0] = gradient wrt x coordinate of atom 0
**     grad[1] = gradient wrt y coordinate of atom 0
**     grad[8] = gradient wrt z coordinate of atom 2
**
**   returns: none
** \ingroup (CHKPT)
*/

void chkpt_wt_grad(double *grad)
{
  int natom;
  char *keyword;
  keyword = chkpt_build_keyword("Energy Gradient");

  natom = chkpt_rd_natom();

  psio_write_entry(PSIF_CHKPT, keyword, (char *) grad, natom*3*sizeof(double));

  free(keyword);
}
