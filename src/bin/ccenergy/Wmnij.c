#include <stdio.h>
#include <libciomr.h>
#include <dpd.h>
#include <qt.h>
#define EXTERN
#include "globals.h"

void Wmnij_build(void)
{
  dpdbuf4 A_anti, A;
  dpdbuf4 WMNIJ, Wmnij, WMnIj, W;
  dpdfile2 tIA, tia;
  dpdbuf4 Eijka, Eijka_anti, Eaijk, Eaijk_anti;
  dpdbuf4 D_anti, D, tauIJAB, tauijab, tauIjAb;

  timer_on("Wmnij");

  dpd_buf4_init(&A_anti, CC_AINTS, 0, 2, 2, 0, 0, 1, "A <ij|kl>");
  dpd_buf4_copy(&A_anti, CC_HBAR, "WMNIJ");
  dpd_buf4_copy(&A_anti, CC_HBAR, "Wmnij");
  dpd_buf4_close(&A_anti);

  dpd_buf4_init(&A, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
  dpd_buf4_copy(&A, CC_HBAR, "WMnIj");
  dpd_buf4_close(&A);
  
  dpd_buf4_init(&WMNIJ, CC_HBAR, 0, 2, 0, 2, 2, 0, "WMNIJ");
  dpd_buf4_init(&Wmnij, CC_HBAR, 0, 2, 0, 2, 2, 0, "Wmnij");
  dpd_buf4_init(&WMnIj, CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");

  dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

  dpd_buf4_init(&Eijka_anti, CC_EINTS, 0, 2, 10, 2, 10, 0,
		"E <ij||ka> (i>j,ka)");
  dpd_buf4_init(&Eijka, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  dpd_buf4_init(&Eaijk_anti, CC_EINTS, 0, 11, 2, 11, 0, 1, "E <ai|jk>");
  dpd_buf4_init(&Eaijk, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");

  dpd_buf4_init(&W, CC_TMP0, 0, 2, 0, 2, 0, 0, "W (MN,IJ)");
  dpd_contract424(&Eijka_anti, &tIA, &W, 3, 1, 0, 1, 0);
  dpd_contract244(&tIA, &Eaijk_anti, &W, 1, 0, 1, 1, 1);
  dpd_buf4_axpy(&W, &WMNIJ, 1);
  dpd_buf4_close(&W);

  dpd_buf4_init(&W, CC_TMP0, 0, 2, 0, 2, 0, 0, "W (MN,IJ)");
  dpd_contract424(&Eijka_anti, &tia, &W, 3, 1, 0, 1, 0);
  dpd_contract244(&tia, &Eaijk_anti, &W, 1, 0, 1, 1, 1);
  dpd_buf4_axpy(&W, &Wmnij, 1);
  dpd_buf4_close(&W);

  dpd_contract424(&Eijka, &tia, &WMnIj, 3, 1, 0, 1, 1);
  dpd_contract244(&tIA, &Eaijk, &WMnIj, 1, 0, 1, 1, 1);

  dpd_buf4_close(&Eijka_anti);
  dpd_buf4_close(&Eijka);
  dpd_buf4_close(&Eaijk_anti);
  dpd_buf4_close(&Eaijk);

  dpd_file2_close(&tIA);
  dpd_file2_close(&tia);

  dpd_buf4_close(&WMNIJ);
  dpd_buf4_close(&Wmnij);
  dpd_buf4_close(&WMnIj);

  dpd_buf4_init(&WMNIJ, CC_HBAR, 0, 2, 2, 2, 2, 0, "WMNIJ");
  dpd_buf4_init(&Wmnij, CC_HBAR, 0, 2, 2, 2, 2, 0, "Wmnij");
  dpd_buf4_init(&WMnIj, CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");

  dpd_buf4_init(&D_anti, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");

  dpd_buf4_init(&tauIJAB, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
  dpd_buf4_init(&tauijab, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
  dpd_buf4_init(&tauIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");

  dpd_contract444(&D_anti, &tauIJAB, &WMNIJ, 0, 0, 1, 1);
  dpd_contract444(&D_anti, &tauijab, &Wmnij, 0, 0, 1, 1);
  dpd_contract444(&D, &tauIjAb, &WMnIj, 0, 0, 1, 1);

  dpd_buf4_close(&tauIJAB);
  dpd_buf4_close(&tauijab);
  dpd_buf4_close(&tauIjAb);

  dpd_buf4_close(&D_anti);
  dpd_buf4_close(&D);

  dpd_buf4_close(&WMNIJ);
  dpd_buf4_close(&Wmnij);
  dpd_buf4_close(&WMnIj);

  timer_off("Wmnij");
}
