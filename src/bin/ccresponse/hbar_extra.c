#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void hbar_extra(void) {
  dpdfile2 lt;
  dpdbuf4 W1, W2, W;
  dpdbuf4 t2, l2;

  /* LIjAb * TIjAb */
  dpd_file2_init(&lt, CC_TMP0, 0, 0, 0, "Lt_IJ");

  dpd_buf4_init(&l2, CC_LAMPS, 0, 0, 7, 2, 7, 0, "LIJAB 0 -1");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
  dpd_contract442(&l2, &t2, &lt, 0, 0, -1.0, 0.0);
  dpd_buf4_close(&t2);
  dpd_buf4_close(&l2);

  dpd_buf4_init(&l2, CC_LAMPS, 0, 0, 5, 0, 5, 0, "LIjAb 0 -1");
  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract442(&l2, &t2, &lt, 0, 0, -1.0, 1.0);
  dpd_buf4_close(&t2);
  dpd_buf4_close(&l2);

  dpd_file2_close(&lt);

  /* 2 W(ME,jb) + W(Me,Jb) */
  dpd_buf4_init(&W1, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
  dpd_buf4_copy(&W1, CC_HBAR, "2 W(ME,jb) + W(Me,Jb)");
  dpd_buf4_close(&W1);
  dpd_buf4_init(&W1, CC_HBAR, 0, 10, 10, 10, 10, 0, "2 W(ME,jb) + W(Me,Jb)");
  dpd_buf4_init(&W2, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
  dpd_buf4_axpy(&W2, &W1, 2);
  dpd_buf4_close(&W2);
  dpd_buf4_sort(&W1, CC_HBAR, rspq, 10, 10, "2 W(jb,ME) + W(Jb,Me)");
  dpd_buf4_close(&W1);

  /* used in WamefSD */
  dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
  dpd_buf4_scmcopy(&W, CC_HBAR, "WAmEf 2(Am,Ef) - (Am,fE)", 2);
  dpd_buf4_sort_axpy(&W, CC_HBAR, pqsr, 11, 5, "WAmEf 2(Am,Ef) - (Am,fE)", -1);
  dpd_buf4_close(&W);

  /*
  dpd_buf4_init(&W1, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe (Mn,eI)");
  dpd_buf4_sort(&W1, CC_HBAR, qprs, 0, 11, "WMnIe - 2WnMIe");
  dpd_buf4_close(&W1);
  dpd_buf4_init(&W1, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe - 2WnMIe");
  dpd_buf4_scm(&W1, -2.0);
  dpd_buf4_init(&W2, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe (Mn,eI)");
  dpd_buf4_axpy(&W2, &W1, 1.0);
  dpd_buf4_close(&W2);
  dpd_buf4_close(&W1);
  dpd_buf4_init(&W1, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe (Mn,eI)");
  dpd_buf4_sort(&W1, CC_HBAR, qprs, 0, 11, "2WMnIe - WnMIe");
  dpd_buf4_close(&W1);
  dpd_buf4_init(&W1, CC_HBAR, 0, 0, 11, 0, 11, 0, "2WMnIe - WnMIe");
  dpd_buf4_scm(&W1, -1.0);
  dpd_buf4_init(&W2, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe (Mn,eI)");
  dpd_buf4_axpy(&W2, &W1, 2.0);
  dpd_buf4_close(&W2);
  dpd_buf4_close(&W1);
  */
}
