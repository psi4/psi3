#include<stdio.h>
#include<stdlib.h>
#include<libint.h>
#ifdef INCLUDE_Default_Deriv1
 #include<libderiv.h>
#endif
#ifdef INCLUDE_MP2R12
 #include<libr12.h>
#endif

#include"defines.h"
#define EXTERN
#include"global.h"

#include"small_fns.h"

/*---------------------------------------------------
  This function checks if requested angular momentum
  is too high for the linked integrals libraries
 ---------------------------------------------------*/
void check_max_am()
{
  if (BasisSet.max_am > LIBINT_MAX_AM || CINTS_MAX_AM < LIBINT_MAX_AM) {
    punt("Angular momentum limit exceeded, link CINTS against a LIBINT library with higher NEW_AM");
  }
  
#ifdef INCLUDE_Default_Deriv1
  if (UserOptions.make_deriv1 && (BasisSet.max_am > LIBDERIV_MAX_AM ||
				  CINTS_MAX_AM < LIBDERIV_MAX_AM + DERIV_LVL) ) {
    punt("Angular momentum limit exceeded, link CINTS against a LIBDERIV library with higher NEW_AM");
  }
#endif

#ifdef INCLUDE_Default_Deriv2
  if (UserOptions.make_deriv2 && (BasisSet.max_am > LIBDERIV_MAX_AM ||
				  CINTS_MAX_AM < LIBDERIV_MAX_AM + DERIV_LVL) ) {
    punt("Angular momentum limit exceeded, link CINTS against a LIBDERIV library with higher NEW_AM");
  }
#endif

#ifdef INCLUDE_MP2R12
  if ((UserOptions.make_r12ints || UserOptions.make_mp2r12) &&
      (BasisSet.max_am > LIBR12_MAX_AM ||
       CINTS_MAX_AM < LIBR12_MAX_AM + 1)) {
    punt("Angular momentum limit exceeded, link CINTS against a LIBR12 library with higher NEW_AM");
  }
#endif

  return;
}
