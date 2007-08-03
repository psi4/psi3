/*! \file 
    \ingroup (CCDENSITY)
    \brief Enter brief description of file here 
*/
#define EXTERN
#include "globals.h"

/* build_Z():  Solve the orbital Z-vector equations:
**
**    sum E,M A(AI,EM) D(orb)(E,M) = -X(A,I)
**
** where A(AI,EM) is the orbital Hessian computed in build_A(), X(A,I)
** is the orbital rotation gradient computed in build_X(), and
** D(orb)(E,M) is the final Z-vector we want. 
**
*/

void build_Z_ROHF(void);
void build_Z_UHF(void);

void build_Z(void)
{
  if(params.ref == 0 || params.ref == 1) build_Z_ROHF();
  else if(params.ref == 2) build_Z_UHF();
}


