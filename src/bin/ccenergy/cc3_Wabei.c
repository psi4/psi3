#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void cc3_Wabei(void)
{
  dpd_buf4_init(&F, CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
  dpd_buf4_copy(&F, CC_MISC, "CC3 WAbEi (Ei,Ab)");
  dpd_buf4_close(&F);


}
