#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void sort_lamps(void)
{
  dpdbuf4 L;

  /* RAK fixing this for new cclambda, assuming A1 ground lambda? */
  dpd_buf4_init(&L, CC_LAMPS, 0, 0, 5, 0, 5, 0, "LIjAb 0 -1");
  dpd_buf4_scmcopy(&L, CC_LAMPS, "2 LIjAb - LIjBa", 2);
  dpd_buf4_sort_axpy(&L, CC_LAMPS, pqsr, 0, 5, "2 LIjAb - LIjBa", -1);
  dpd_buf4_close(&L);

  /*
  dpd_buf4_init(&L, CC_LAMPS, 0, 0, 5, 0, 5, 0, "LIjAb");
  dpd_buf4_scmcopy(&L, CC_LAMPS, "2 LIjAb - LIjBa", 2);
  dpd_buf4_sort_axpy(&L, CC_LAMPS, pqsr, 0, 5, "2 LIjAb - LIjBa", -1);
  dpd_buf4_close(&L);
  */
}  
