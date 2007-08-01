/*! \file oeprop_ints.c
    \ingroup (CINTS)
    \brief Enter brief description of file here 
*/
#include<stdio.h>
#include<libint/libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"

#include"moment_ints.h"
#include"moment_deriv1.h"
#include"angmom_ints.h"

void oeprop_ints()
{
  moment_ints();
  moment_deriv1();
  angmom_ints();
  return;
}
