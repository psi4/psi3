#include <stdio.h>
#include <libciomr.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void Wmnij_build(void)
{
  struct dpdbuf A_anti, A;
  struct dpdbuf WMNIJ, Wmnij, WMnIj;
  struct oe_dpdfile tIA, tia;
  struct dpdbuf Eijka, Eijka_anti, Eaijk, Eaijk_anti;
  struct dpdbuf D_anti, D, tauIJAB, tauijab, tauIjAb;

  dpd_buf_init(&A_anti, CC_AINTS, 2, 2, 0, 0, 1, "A <ij|kl>", 0, outfile);
  dpd_copy(&A_anti, CC_HBAR, "WMNIJ", 0, outfile);
  dpd_copy(&A_anti, CC_HBAR, "Wmnij", 0, outfile);
  dpd_buf_close(&A_anti);

  dpd_buf_init(&A, CC_AINTS, 0, 0, 0, 0, 0, "A <ij|kl>", 0, outfile);
  dpd_copy(&A, CC_HBAR, "WMnIj", 0, outfile);
  dpd_buf_close(&A);
  
  dpd_buf_init(&WMNIJ, CC_HBAR, 2, 0, 2, 2, 0, "WMNIJ", 0, outfile);
  dpd_buf_init(&Wmnij, CC_HBAR, 2, 0, 2, 2, 0, "Wmnij", 0, outfile);
  dpd_buf_init(&WMnIj, CC_HBAR, 0, 0, 0, 0, 0, "WMnIj", 0, outfile);

  dpd_oe_file_init(&tIA, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_init(&tia, CC_OEI, 0, 1, "tia", 0, outfile);

  dpd_buf_init(&Eijka_anti, CC_EINTS, 2, 10, 2, 10, 0,
	       "E <ij||ka> (i>j,ka)", 0, outfile);
  dpd_buf_init(&Eijka, CC_EINTS, 0, 10, 0, 10, 0, "E <ij|ka>", 0, outfile);
  dpd_buf_init(&Eaijk_anti,CC_EINTS, 11, 2, 11, 0, 1, "E <ai|jk>", 0, outfile);
  dpd_buf_init(&Eaijk, CC_EINTS, 11, 0, 11, 0, 0, "E <ai|jk>", 0, outfile);

  dpd_contract221(&Eijka_anti, &tIA, &WMNIJ, 3, 1, 0, 1, 1, 0, outfile);
  dpd_contract221(&Eijka_anti, &tia, &Wmnij, 3, 1, 0, 1, 1, 0, outfile);
  dpd_contract221(&Eijka, &tia, &WMnIj, 3, 1, 0, 1, 1, 0, outfile);

  dpd_contract212(&tIA, &Eaijk_anti, &WMNIJ, 1, 0, 1, 1, 1, 0, outfile);
  dpd_contract212(&tia, &Eaijk_anti, &Wmnij, 1, 0, 1, 1, 1, 0, outfile);
  dpd_contract212(&tIA, &Eaijk, &WMnIj, 1, 0, 1, 1, 1, 0, outfile);

  dpd_buf_close(&Eijka_anti);
  dpd_buf_close(&Eijka);
  dpd_buf_close(&Eaijk_anti);
  dpd_buf_close(&Eaijk);

  dpd_oe_file_close(&tIA);
  dpd_oe_file_close(&tia);

  dpd_buf_close(&WMNIJ);
  dpd_buf_close(&Wmnij);
  dpd_buf_close(&WMnIj);

  dpd_buf_init(&WMNIJ, CC_HBAR, 2, 2, 2, 2, 0, "WMNIJ", 0, outfile);
  dpd_buf_init(&Wmnij, CC_HBAR, 2, 2, 2, 2, 0, "Wmnij", 0, outfile);
  dpd_buf_init(&WMnIj, CC_HBAR, 0, 0, 0, 0, 0, "WMnIj", 0, outfile);

  dpd_buf_init(&D_anti, CC_DINTS, 2, 7, 2, 7, 0,
	       "D <ij||ab> (i>j,a>b)", 0, outfile);
  dpd_buf_init(&D, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);

  dpd_buf_init(&tauIJAB, CC_TAMPS, 2, 7, 2, 7, 0, "tauIJAB", 0, outfile);
  dpd_buf_init(&tauijab, CC_TAMPS, 2, 7, 2, 7, 0, "tauijab", 0, outfile);
  dpd_buf_init(&tauIjAb, CC_TAMPS, 0, 5, 0, 5, 0, "tauIjAb", 0, outfile);

  dpd_contract222(&D_anti, &tauIJAB, &WMNIJ, 0, 0, 1, 1, 0, outfile);
  dpd_contract222(&D_anti, &tauijab, &Wmnij, 0, 0, 1, 1, 0, outfile);
  dpd_contract222(&D, &tauIjAb, &WMnIj, 0, 0, 1, 1, 0, outfile);

  dpd_buf_close(&tauIJAB);
  dpd_buf_close(&tauijab);
  dpd_buf_close(&tauIjAb);

  dpd_buf_close(&D_anti);
  dpd_buf_close(&D);

  dpd_buf_close(&WMNIJ);
  dpd_buf_close(&Wmnij);
  dpd_buf_close(&WMnIj);
}
