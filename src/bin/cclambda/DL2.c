#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void DL2(void)
{
  dpdbuf4 D;

  /* RHS = <ij||ab> */
  dpd_buf4_init(&D, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
  dpd_buf4_copy(&D, CC_LAMPS, "New LIJAB");
  dpd_buf4_copy(&D, CC_LAMPS, "New Lijab");
  dpd_buf4_close(&D);

  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_copy(&D, CC_LAMPS, "New LIjAb");
  dpd_buf4_close(&D);
}
