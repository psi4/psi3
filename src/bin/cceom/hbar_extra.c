#include <stdio.h>
#include <math.h>
#define EXTERN
#include "globals.h"

void hbar_extra(void) {
  dpdbuf4 W, W1, W2, WAmEf, WmBeJ, WmBEj, WmNIe, WMnIe;

  if (params.eom_ref == 2) {
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 20, 20, 20, 20, 0, "WMBEJ");
    dpd_buf4_sort(&W, CC_HBAR, rspq, 20, 20, "WMBEJ (JB,ME)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 30, 20, 30, 20, 0, "WmBeJ"); /* (me,JB) */
    dpd_buf4_sort(&W, CC_HBAR, rspq, 20, 30, "WmBeJ (JB,me)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 30, 30, 30, 30, 0, "Wmbej");
    dpd_buf4_sort(&W, CC_HBAR, rspq, 30, 30, "Wmbej (jb,me)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 20, 30, 20, 30, 0, "WMbEj"); /* (ME,jb) */
    dpd_buf4_sort(&W, CC_HBAR, rspq, 30, 20, "WMbEj (jb,ME)");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_HBAR, H_IRR, 27, 23, 27, 23, 0, "WmBiJ");
    dpd_buf4_sort(&W, CC_HBAR, pqsr, 27, 22, "WmBiJ (mB,Ji)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 27, 22, 27, 22, 0, "WmBiJ (mB,Ji)");
    dpd_buf4_sort(&W, CC_HBAR, qprs, 26, 22, "WmBiJ (Bm,Ji)");
    dpd_buf4_close(&W);
  }

  if ((params.eom_ref == 0) || (params.eom_ref == 1)) {
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
  }

  if (params.eom_ref == 0 ) { /* RHF */
    /* 2 W(ME,jb) + W(Me,Jb) */
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 10, 10, 10, 10, 0, "WMbeJ");
    dpd_buf4_copy(&W, CC_HBAR, "2 W(ME,jb) + W(Me,Jb)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W1, CC_HBAR, H_IRR, 10, 10, 10, 10, 0, "2 W(ME,jb) + W(Me,Jb)");
    dpd_buf4_init(&W2, CC_HBAR, H_IRR, 10, 10, 10, 10, 0, "WMbEj");
    dpd_buf4_axpy(&W2, &W1, 2);
    dpd_buf4_close(&W2);
    dpd_buf4_close(&W1);

    /* (2WmBeJ + WmAEi) (jb,me) */
    dpd_buf4_init(&WmBeJ, CC_HBAR, H_IRR, 10, 10, 10, 10, 0, "WmBeJ"); /* (me,JB) */
    dpd_buf4_init(&WmBEj, CC_HBAR, H_IRR, 10, 10, 10, 10, 0, "WmBEj"); /* (ME,jb) */
    dpd_buf4_copy(&WmBeJ, CC_HBAR, "(2WmBeJ + WmBEj) (me,jb)");
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 10, 10, 10, 10, 0, "(2WmBeJ + WmBEj) (me,jb)");
    dpd_buf4_axpy(&WmBeJ, &W, 1.0);
    dpd_buf4_axpy(&WmBEj, &W, 1.0);
    dpd_buf4_close(&WmBeJ);
    dpd_buf4_close(&WmBEj);
    dpd_buf4_sort(&W, CC_HBAR, rspq, 10, 10, "(2WmBeJ + WmBEj) (jb,me)");
    dpd_buf4_close(&W);

    /* used in RHF WmnieSD */
    dpd_buf4_init(&WMnIe, CC_HBAR, H_IRR, 0, 11, 0, 11, 0, "WMnIe");
    dpd_buf4_sort(&WMnIe, CC_HBAR, qprs, 0, 11, "WMnIe - 2WnMIe");
    dpd_buf4_close(&WMnIe);
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 0, 11, 0, 11, 0, "WMnIe - 2WnMIe");
    dpd_buf4_scm(&W, -2.0);
    dpd_buf4_init(&WMnIe, CC_HBAR, H_IRR, 0, 11, 0, 11, 0, "WMnIe");
    dpd_buf4_axpy(&WMnIe, &W, 1.0);
    dpd_buf4_close(&WMnIe);
    dpd_buf4_close(&W);

    /* used in RHF WnmjeDS */
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 0, 11, 0, 11, 0, "WMnIe");
    dpd_buf4_sort(&W, CC_HBAR, qprs, 0, 11, "2WMnIe - WnMIe");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 0, 11, 0, 11, 0, "2WMnIe - WnMIe");
    dpd_buf4_scm(&W, -1.0);
    dpd_buf4_init(&WMnIe, CC_HBAR, H_IRR, 0, 11, 0, 11, 0, "WMnIe");
    dpd_buf4_axpy(&WMnIe, &W, 2.0);
    dpd_buf4_close(&WMnIe);
    dpd_buf4_close(&W);

    /* used in WamefSD */
    dpd_buf4_init(&WAmEf, CC_HBAR, H_IRR, 11, 5, 11, 5, 0, "WAmEf (Am,Ef)");
    dpd_buf4_sort(&WAmEf, CC_HBAR, pqsr, 11, 5, "WAmEf 2(Am,Ef) - (Am,fE)");
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)");
    dpd_buf4_scm(&W, -1.0);
    dpd_buf4_axpy(&WAmEf, &W, 2.0); 
    dpd_buf4_close(&WAmEf);
    dpd_buf4_close(&W);

    /* used in WbmfeDS */
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)");
    dpd_buf4_sort(&W, CC_HBAR, qprs, 10, 5, "WAmEf 2(mA,Ef) - (mA,fE)");
    dpd_buf4_close(&W);

  }

  return;
}
