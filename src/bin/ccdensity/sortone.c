/*! \file 
    \ingroup (CCDENSITY)
    \brief Enter brief description of file here 
*/
#define EXTERN
#include "globals.h"

void sortone(struct RHO_Params rho_params)
{
  if(params.ref == 0 || params.ref == 1) sortone_ROHF(rho_params);
  else if(params.ref == 2) sortone_UHF(rho_params);
}
