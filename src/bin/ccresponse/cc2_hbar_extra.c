#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void cc2_hbar_extra(void) {
  dpdfile2 t1;
  dpdbuf4 A, D, E, Z, Z1;
  dpdbuf4 W1, W2, W;

  /* 2 W(ME,jb) + W(Me,Jb) */
  dpd_buf4_init(&W1, CC_HBAR, 0, 10, 10, 10, 10, 0, "CC2 WMbeJ");
  dpd_buf4_copy(&W1, CC_HBAR, "CC2 2 W(ME,jb) + W(Me,Jb)");
  dpd_buf4_close(&W1);
  dpd_buf4_init(&W1, CC_HBAR, 0, 10, 10, 10, 10, 0, "CC2 2 W(ME,jb) + W(Me,Jb)");
  dpd_buf4_init(&W2, CC_HBAR, 0, 10, 10, 10, 10, 0, "CC2 WMbEj");
  dpd_buf4_axpy(&W2, &W1, 2);
  dpd_buf4_close(&W2);
  dpd_buf4_sort(&W1, CC_HBAR, rspq, 10, 10, "CC2 2 W(jb,ME) + W(Jb,Me)");
  dpd_buf4_close(&W1);

  /* used in WamefSD */
  dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
  dpd_buf4_scmcopy(&W, CC_HBAR, "WAmEf 2(Am,Ef) - (Am,fE)", 2);
  dpd_buf4_sort_axpy(&W, CC_HBAR, pqsr, 11, 5, "WAmEf 2(Am,Ef) - (Am,fE)", -1);
  dpd_buf4_close(&W);

  /* CC2 WMnIj */

  dpd_buf4_init(&A, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
  dpd_buf4_copy(&A, CC_HBAR, "CC2 WMnIj (Mn,Ij)");
  dpd_buf4_close(&A);

  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");

  /* Wmnij <- + P(ij) t(j,e) * <mn||ie> */
  dpd_buf4_init(&Z, CC_TMP0, 0, 0, 0, 0, 0, 0, "CC2 ZMnIj (Mn,Ij)");
  dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  dpd_contract424(&E, &t1, &Z, 3, 1, 0, 1, 0);
  dpd_buf4_close(&E);

  dpd_buf4_init(&W, CC_HBAR, 0, 0, 0, 0, 0, 0, "CC2 WMnIj (Mn,Ij)");
  dpd_buf4_axpy(&Z, &W, 1);
  dpd_buf4_close(&W);
  dpd_buf4_sort_axpy(&Z, CC_HBAR, qpsr, 0, 0, "CC2 WMnIj (Mn,Ij)", 1);
  dpd_buf4_close(&Z);

  /* Wmnij<- +1/2 P(ij) t(i,e) t(j,f) * <mn||ef> */
  dpd_buf4_init(&Z, CC_TMP0, 0, 0, 10, 0, 10, 0, "CC2 ZMnIf (Mn,If)");
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_contract244(&t1, &D, &Z, 1, 2, 1, 1, 0);
  dpd_buf4_close(&D);

  dpd_buf4_init(&Z1, CC_TMP0, 0, 0, 0, 0, 0, 0, "CC2 ZMnIj (Mn,Ij)");
  dpd_contract424(&Z, &t1, &Z1, 3, 1, 0, 0.5, 0);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&W, CC_HBAR, 0, 0, 0, 0, 0, 0, "CC2 WMnIj (Mn,Ij)");
  dpd_buf4_axpy(&Z1, &W, 1);
  dpd_buf4_close(&W);
  dpd_buf4_sort_axpy(&Z1, CC_HBAR, qpsr, 0, 0, "CC2 WMnIj (Mn,Ij)", 1);
  dpd_buf4_close(&Z1);

  dpd_file2_close(&t1);
}
