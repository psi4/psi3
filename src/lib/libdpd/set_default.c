/*! \file 
    \ingroup (DPD)
    \brief Enter brief description of file here 
*/
#include "dpd.h"
#define EXTERN
#include "dpd.gbl"

int dpd_set_default(int dpd_num)
{
  dpd_default = dpd_num;

  return 0;
}
