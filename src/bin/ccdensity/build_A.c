/*! \file build_A.c
    \ingroup (CCDENSITY)
    \brief Enter brief description of file here 
*/
#define EXTERN
#include "globals.h"

/* BUILD_A(): Construct the molecular orbital Hessian, A.
** */

void build_A_ROHF(void);
void build_A_UHF(void);

void build_A(void)
{
  if(params.ref == 0 || params.ref == 1) build_A_ROHF();
  else if(params.ref == 2) build_A_UHF();
}

