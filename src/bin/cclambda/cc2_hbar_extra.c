#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void cc2_hbar_extra(void) {
  dpdbuf4 W1, W2, W;

  if(params.ref == 0) { /** RHF **/
    /* 2 W(ME,jb) + W(Me,Jb) */
    dpd_buf4_init(&W1, CC_HBAR, 0, 10, 10, 10, 10, 0, "CC2 WMbeJ");
    dpd_buf4_copy(&W1, CC_HBAR, "CC2 2 W(ME,jb) + W(Me,Jb)");
    dpd_buf4_close(&W1);
    dpd_buf4_init(&W1, CC_HBAR, 0, 10, 10, 10, 10, 0, "CC2 2 W(ME,jb) + W(Me,Jb)");
    dpd_buf4_init(&W2, CC_HBAR, 0, 10, 10, 10, 10, 0, "CC2 WMbEj");
    dpd_buf4_axpy(&W2, &W1, 2);
    dpd_buf4_close(&W2);
    dpd_buf4_close(&W1);

    dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "CC2 WAmEf");
    dpd_buf4_scmcopy(&W, CC_HBAR, "CC2 WAmEf 2(Am,Ef) - (Am,fE)", 2);
    dpd_buf4_sort_axpy(&W, CC_HBAR, pqsr, 11, 5, "CC2 WAmEf 2(Am,Ef) - (Am,fE)", -1);
    dpd_buf4_close(&W);
  }
}
