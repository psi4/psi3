#include "dpd.h"
#define EXTERN
#include "dpd.gbl"

int dpd_set_default(int dpd_num)
{
  dpd_default = &(dpd_list[dpd_num]);

  return 0;
}
