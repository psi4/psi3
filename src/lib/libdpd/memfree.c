/*
** Function to return number of double words available for allocation.
*/

#include <stdio.h>
#include "dpd.h"
#define EXTERN
#include "dpd.gbl"

int dpd_memfree(void)
{
  return dpd_default->memory - (dpd_default->memused - 
				dpd_default->memcache + 
				dpd_default->memlocked);
}
