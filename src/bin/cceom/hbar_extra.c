#include <stdio.h>
#include <math.h>
#define EXTERN
#include "globals.h"

void hbar_extra(void) {
  dpdbuf4 W;

  dpd_buf4_init(&W, CC_HBAR, H_IRR, 10, 10, 10, 10, 0, "WMBEJ");
  dpd_buf4_sort(&W, CC_HBAR, rspq, 10, 10, "WMBEJ (JB,ME)");
  dpd_buf4_close(&W);
  dpd_buf4_init(&W, CC_HBAR, H_IRR, 10, 10, 10, 10, 0, "WmBeJ");
  dpd_buf4_sort(&W, CC_HBAR, rspq, 10, 10, "WmBeJ (JB,me)");
  dpd_buf4_close(&W);
  dpd_buf4_init(&W, CC_HBAR, H_IRR, 10, 10, 10, 10, 0, "Wmbej");
  dpd_buf4_sort(&W, CC_HBAR, rspq, 10, 10, "Wmbej (jb,me)");
  dpd_buf4_close(&W);
  dpd_buf4_init(&W, CC_HBAR, H_IRR, 10, 10, 10, 10, 0, "WMbEj");
  dpd_buf4_sort(&W, CC_HBAR, rspq, 10, 10, "WMbEj (jb,ME)");
  dpd_buf4_close(&W);

  dpd_buf4_init(&W, CC_HBAR, H_IRR, 10, 7, 10, 7, 0, "WAMEF");
  dpd_buf4_sort(&W, CC_HBAR, qprs, 11, 7, "WAMEF (AM,E>F)");
  dpd_buf4_close(&W);
  dpd_buf4_init(&W, CC_HBAR, H_IRR, 10, 7, 10, 7, 0, "Wamef");
  dpd_buf4_sort(&W, CC_HBAR, qprs, 11, 7, "Wamef (am,e>f)");
  dpd_buf4_close(&W);
  dpd_buf4_init(&W, CC_HBAR, H_IRR, 10, 5, 10, 5, 0, "WAmEf");
  dpd_buf4_sort(&W, CC_HBAR, qprs, 11, 5, "WAmEf (Am,Ef)");
  dpd_buf4_close(&W);
  dpd_buf4_init(&W, CC_HBAR, H_IRR, 10, 5, 10, 5, 0, "WaMeF");
  dpd_buf4_sort(&W, CC_HBAR, qprs, 11, 5, "WaMeF (aM,eF)");
  dpd_buf4_close(&W);

  dpd_buf4_init(&W, CC_HBAR, H_IRR, 10, 0, 10, 0, 0, "WmBiJ");
  dpd_buf4_sort(&W, CC_HBAR, pqsr, 10, 0, "WmBiJ (mB,Ji)");
  dpd_buf4_close(&W);
  dpd_buf4_init(&W, CC_HBAR, H_IRR, 10, 0, 10, 0, 0, "WmBiJ (mB,Ji)");
  dpd_buf4_sort(&W, CC_HBAR, qprs, 11, 0, "WmBiJ (Bm,Ji)");
  dpd_buf4_close(&W);

  dpd_buf4_init(&W, CC_HBAR, H_IRR, 11, 5, 11, 5, 0, "WeIaB");
  dpd_buf4_sort(&W, CC_HBAR, qprs, 10, 5, "WeIaB (Ie,aB)");
  dpd_buf4_close(&W);
  dpd_buf4_init(&W, CC_HBAR, H_IRR, 10, 5, 10, 5, 0, "WeIaB (Ie,aB)");
  dpd_buf4_sort(&W, CC_HBAR, pqsr, 10, 5, "WeIaB (Ie,Ab)");
  dpd_buf4_close(&W);

  return;
}
