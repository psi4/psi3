#include<stdio.h>
#include<libint/libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"

#include"moment_ints.h"
#include"moment_deriv1.h"
#include"angmom_ints.h"
#include"enm_deriv_ints.h"

void oeprop_ints()
{
  moment_ints();
  moment_deriv1();
  angmom_ints();
  enm_deriv_ints();
  return;
}
