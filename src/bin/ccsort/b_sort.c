#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

void b_sort(void)
{
  dpdbuf4 B;

  /* <ab||cd> */
  if(params.ref != 2) { /* RHF/ROHF references */
    dpd_buf4_init(&B, CC_BINTS, 0, 7, 7, 5, 5, 1, "B <ab|cd>");
    dpd_buf4_copy(&B, CC_BINTS, "B <ab||cd> (a>b,c>d)");
    dpd_buf4_close(&B);
  }
}
