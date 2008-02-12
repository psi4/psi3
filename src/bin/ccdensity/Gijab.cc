/*! \file Gijab.cc
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void Gijab_ROHF(void);
void Gijab_UHF(void);

void Gijab(void)
{
  if(params.ref == 0 || params.ref == 1) Gijab_ROHF();
  else if(params.ref == 2) Gijab_UHF();
}

}} // namespace psi::ccdensity
