#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void DT2(void)
{
  dpdbuf4 D;

  if(params.ref == 0) { /*** RHF ***/
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_buf4_copy(&D, CC_TAMPS, "New tIjAb");
    dpd_buf4_close(&D);
  }
  else if(params.ref == 1) { /*** ROHF ***/

    dpd_buf4_init(&D, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
    dpd_buf4_copy(&D, CC_TAMPS, "New tIJAB");
    dpd_buf4_copy(&D, CC_TAMPS, "New tijab");
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_buf4_copy(&D, CC_TAMPS, "New tIjAb");
    dpd_buf4_close(&D);
  }
}
