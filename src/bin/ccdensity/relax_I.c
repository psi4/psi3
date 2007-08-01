/*! \file relax_I.c
    \ingroup (CCDENSITY)
    \brief Enter brief description of file here 
*/
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/* relax_I(): Add the orbital-response contributions from the
** one-electron density matrix to the I(I,J) and I(I,A) blocks of the
** Lagrangian.  These terms arise from the first-order CPHF
** equations.  I *think* the following code is general enough to deal
** with both RHF and ROHF cases.
** */

void relax_I_ROHF(void);
void relax_I_UHF(void);

void relax_I(void)
{
  if(params.ref == 0 || params.ref == 1) relax_I_ROHF();
  else if(params.ref == 2) relax_I_UHF();
}
 
