#include <stdio.h>
#include <math.h>
#define EXTERN
#include "globals.h"

/* This function computes the H-bar doubles-singles block contribution
of Wabej to a Sigma vector stored at Sigma plus 'i' */

void WabejDS(int i, int irrep) {
  dpdfile2 CME, Cme;
  dpdbuf4 SIJAB, Sijab, SIjAb;
  dpdbuf4 WABEI, Wabei, WAbEi, WaBeI, WM, WP;
  char CME_lbl[32], Cme_lbl[32], SIJAB_lbl[32], Sijab_lbl[32], SIjAb_lbl[32];

  sprintf(CME_lbl, "%s %d", "CME", i);
  sprintf(Cme_lbl, "%s %d", "Cme", i);
  sprintf(SIJAB_lbl, "%s %d", "SIJAB", i);
  sprintf(Sijab_lbl, "%s %d", "Sijab", i);
  sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

  /* SIJAB += WABEJ * CIE - WABEI * CJE */
  dpd_buf4_init(&WP, EOM_TMP, irrep, 0, 7, 0, 7, 0, "WabejDS_P");
  dpd_buf4_init(&WABEI, CC_HBAR, irrep, 11, 7, 11, 7, 0, "WEIAB");
  dpd_file2_init(&CME, EOM_CME, irrep, 0, 1, CME_lbl);
  dpd_contract244(&CME, &WABEI, &WP, 1, 0, 0, 1.0, 0.0);
  dpd_file2_close(&CME);
  dpd_buf4_close(&WABEI);
  dpd_buf4_sort(&WP, EOM_TMP, qprs, 0, 7, "WabejDS_M");
  dpd_buf4_init(&SIJAB, EOM_SIJAB, irrep, 0, 7, 2, 7, 0, SIJAB_lbl);
  dpd_buf4_axpy(&WP, &SIJAB, 1.0);
  dpd_buf4_close(&WP);
  dpd_buf4_init(&WM, EOM_TMP, irrep, 0, 7, 0, 7, 0, "WabejDS_M");
  dpd_buf4_axpy(&WM, &SIJAB, -1.0);
  dpd_buf4_close(&WM);
  dpd_buf4_close(&SIJAB);

  /* Sijab += Wabej * Cie - Wabei * Cje */
  dpd_buf4_init(&WP, EOM_TMP, irrep, 0, 7, 0, 7, 0, "WabejDS_P");
  dpd_buf4_init(&Wabei, CC_HBAR, irrep, 11, 7, 11, 7, 0, "Weiab");
  dpd_file2_init(&Cme, EOM_Cme, irrep, 0, 1, Cme_lbl);
  dpd_contract244(&Cme, &Wabei, &WP, 1, 0, 0, 1.0, 0.0);
  dpd_file2_close(&Cme);
  dpd_buf4_close(&Wabei);
  dpd_buf4_sort(&WP, EOM_TMP, qprs, 0, 7, "WabejDS_M");
  dpd_buf4_init(&Sijab, EOM_Sijab, irrep, 0, 7, 2, 7, 0, Sijab_lbl);
  dpd_buf4_axpy(&WP, &Sijab, 1.0);
  dpd_buf4_close(&WP);
  dpd_buf4_init(&WM, EOM_TMP, irrep, 0, 7, 0, 7, 0, "WabejDS_M");
  dpd_buf4_axpy(&WM, &Sijab, -1.0);
  dpd_buf4_close(&WM);
  dpd_buf4_close(&Sijab);


  /* SIjAb += WAbEj * CIE - WAbeI * Cje */
  dpd_buf4_init(&SIjAb, EOM_SIjAb, irrep, 0, 5, 0, 5, 0, SIjAb_lbl);
  dpd_file2_init(&CME, EOM_CME, irrep, 0, 1, CME_lbl);
  dpd_buf4_init(&WAbEi, CC_HBAR, irrep, 11, 5, 11, 5, 0, "WEiAb");
  dpd_contract244(&CME, &WAbEi, &SIjAb, 1, 0, 0, 1.0, 1.0);
  dpd_file2_close(&CME);
  dpd_buf4_close(&WAbEi);
  dpd_buf4_init(&WaBeI, CC_HBAR, irrep, 10, 5, 10, 5, 0, "WeIaB (Ie,Ab)");
  dpd_file2_init(&Cme, EOM_Cme, irrep, 0, 1, Cme_lbl);
  dpd_contract424(&WaBeI, &Cme, &SIjAb, 1, 1, 1, 1.0, 1.0);
  dpd_file2_close(&Cme);
  dpd_buf4_close(&WaBeI);
  dpd_buf4_close(&SIjAb);

#ifdef EOM_DEBUG
  check_sum("WabejDS",i,irrep);
#endif

  return;
}

