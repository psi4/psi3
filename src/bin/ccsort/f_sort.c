#include <dpd.h>
#include "Params.h"
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

void f_sort(void)
{
  dpdbuf4 F;

  /* <ia||bc> */
  dpd_buf4_init(&F, CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
  dpd_buf4_copy(&F, CC_FINTS, "F <ia||bc> (ia,b>c)");
  dpd_buf4_close(&F);

  /* <ai|bc> */
  dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_buf4_sort(&F, CC_FINTS, qpsr, 11, 5, "F <ai|bc>");
  dpd_buf4_close(&F);
}
