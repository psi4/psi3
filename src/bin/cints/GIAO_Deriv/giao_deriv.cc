
#include "giao_oe_deriv.h"
#include "giao_te_deriv.h"

extern "C" void giao_deriv()
{
  giao_oe_deriv();
  giao_te_deriv();
}

