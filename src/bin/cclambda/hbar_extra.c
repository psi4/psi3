#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void hbar_extra(void) {
  dpdbuf4 W1, W2;

  /* 2 W(ME,jb) + W(Me,Jb) */
  dpd_buf4_init(&W1, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
  dpd_buf4_copy(&W1, CC_HBAR, "2 W(ME,jb) + W(Me,Jb)");
  dpd_buf4_close(&W1);
  dpd_buf4_init(&W1, CC_HBAR, 0, 10, 10, 10, 10, 0, "2 W(ME,jb) + W(Me,Jb)");
  dpd_buf4_init(&W2, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
  dpd_buf4_axpy(&W2, &W1, 2);
  dpd_buf4_close(&W2);
  dpd_buf4_close(&W1);

  /* (2WmBeJ + WmBEj) (ME,jb) */
  dpd_buf4_init(&W1, CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBeJ"); /* (me,JB) */
  dpd_buf4_copy(&W1, CC_HBAR, "(2WmBeJ + WmBEj) (me,jb)");
  dpd_buf4_close(&W1);

  dpd_buf4_init(&W1, CC_HBAR, 0, 10, 10, 10, 10, 0, "(2WmBeJ + WmBEj) (me,jb)");
  dpd_buf4_scm(&W1, 2.0);
  dpd_buf4_init(&W2, CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBEj"); /* (ME,jb) */
  dpd_buf4_axpy(&W2, &W1, 1.0);
  dpd_buf4_close(&W2);
  dpd_buf4_sort(&W1, CC_HBAR, rspq, 10, 10, "(2WmBeJ + WmBEj) (jb,me)");
  dpd_buf4_close(&W1);

  /* used in WamefSD */
  dpd_buf4_init(&W1, CC_HBAR, 0, 10, 5, 10, 5, 0, "WAmEf");
  dpd_buf4_sort(&W1, CC_HBAR, qprs, 11, 5, "WAmEf (Am,Ef)");
  dpd_buf4_close(&W1);
  dpd_buf4_init(&W1, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf (Am,Ef)");
  dpd_buf4_sort(&W1, CC_HBAR, pqsr, 11, 5, "WAmEf 2(Am,Ef) - (Am,fE)");
  dpd_buf4_init(&W2, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)");
  dpd_buf4_scm(&W2, -1.0);
  dpd_buf4_axpy(&W1, &W2, 2.0); 
  dpd_buf4_close(&W1);
  dpd_buf4_close(&W2);

  /* used in WbmfeDS */
  dpd_buf4_init(&W1, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)");
  dpd_buf4_sort(&W1, CC_HBAR, qprs, 10, 5, "WAmEf 2(mA,Ef) - (mA,fE)");
  dpd_buf4_close(&W1);

  /* used in RHF WmnieSD */
  dpd_buf4_init(&W1, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe");
  dpd_buf4_sort(&W1, CC_HBAR, qprs, 0, 11, "WMnIe - 2WnMIe");
  dpd_buf4_close(&W1);
  dpd_buf4_init(&W1, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe - 2WnMIe");
  dpd_buf4_scm(&W1, -2.0);
  dpd_buf4_init(&W2, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe");
  dpd_buf4_axpy(&W2, &W1, 1.0);
  dpd_buf4_close(&W2);
  dpd_buf4_close(&W1);

  /* used in RHF WnmjeDS */
  dpd_buf4_init(&W1, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe");
  dpd_buf4_sort(&W1, CC_HBAR, qprs, 0, 11, "2WMnIe - WnMIe");
  dpd_buf4_close(&W1);
  dpd_buf4_init(&W1, CC_HBAR, 0, 0, 11, 0, 11, 0, "2WMnIe - WnMIe");
  dpd_buf4_scm(&W1, -1.0);
  dpd_buf4_init(&W2, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe");
  dpd_buf4_axpy(&W2, &W1, 2.0);
  dpd_buf4_close(&W2);
  dpd_buf4_close(&W1);
}
