#include <dpd.h>
#include "Params.h"
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

void f_sort(void)
{
  dpdbuf4 F;

  if(params.ref == 2) {  /*** UHF ***/
    dpd_buf4_init(&F, CC_FINTS, 0, 28, 26, 28, 26, 0, "F <Ab|Ci>");
    dpd_buf4_sort(&F, CC_FINTS, spqr, 27, 29, "F <iA|bC>");
    dpd_buf4_close(&F);
  }
  else { /*** RHF/ROHF ***/
    /* <ia||bc> */
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
    dpd_buf4_copy(&F, CC_FINTS, "F <ia||bc> (ia,b>c)");
    dpd_buf4_close(&F);

    /* <ai|bc> */
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_buf4_sort(&F, CC_FINTS, qpsr, 11, 5, "F <ai|bc>");
    dpd_buf4_close(&F);
  }
}
