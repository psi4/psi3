/*! \file 
    \ingroup (CCDENSITY)
    \brief Enter brief description of file here 
*/
#define EXTERN
#include "globals.h"

void Gijab_ROHF(void);
void Gijab_UHF(void);

void Gijab(void)
{
  if(params.ref == 0 || params.ref == 1) Gijab_ROHF();
  else if(params.ref == 2) Gijab_UHF();
}
