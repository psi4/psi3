#include <stdio.h>
#include <math.h>
#define EXTERN
#include "globals.h"

/* This function computes the H-bar doubles-singles block contribution
   of Wabej to a Sigma vector stored at Sigma plus 'i' */

void WabejDS(int i, int C_irr) {
  dpdfile2 CME, Cme;
  dpdbuf4 SIJAB, Sijab, SIjAb;
  dpdbuf4 WABEI, Wabei, WAbEi, WaBeI, WM, WP, Z;
  char CME_lbl[32], Cme_lbl[32], SIJAB_lbl[32], Sijab_lbl[32], SIjAb_lbl[32];

  if (params.eom_ref == 0) { /* RHF */
    sprintf(CME_lbl, "%s %d", "CME", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    dpd_buf4_init(&Z, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "WabejDS Z(Ij,Ab)");
    dpd_buf4_init(&WAbEi, CC_HBAR, H_IRR, 11, 5, 11, 5, 0, "WEiAb");
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_contract244(&CME, &WAbEi, &Z, 1, 0, 0, 1.0, 0.0);
    dpd_file2_close(&CME);
    dpd_buf4_close(&WAbEi);
    dpd_buf4_sort(&Z, EOM_TMP, qpsr, 0, 5, "WabejDS Z(jI,bA)");
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    dpd_buf4_axpy(&Z, &SIjAb, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "WabejDS Z(jI,bA)");
    dpd_buf4_axpy(&Z, &SIjAb, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&SIjAb);
  }

  else if (params.eom_ref == 1) { // ROHF
    sprintf(CME_lbl, "%s %d", "CME", i);
    sprintf(Cme_lbl, "%s %d", "Cme", i);
    sprintf(SIJAB_lbl, "%s %d", "SIJAB", i);
    sprintf(Sijab_lbl, "%s %d", "Sijab", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    /* SIJAB += WABEJ * CIE - WABEI * CJE */
    dpd_buf4_init(&WP, EOM_TMP, C_irr, 0, 7, 0, 7, 0, "WabejDS_P");
    dpd_buf4_init(&WABEI, CC_HBAR, H_IRR, 11, 7, 11, 7, 0, "WEIAB");
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_contract244(&CME, &WABEI, &WP, 1, 0, 0, 1.0, 0.0);
    dpd_file2_close(&CME);
    dpd_buf4_close(&WABEI);
    dpd_buf4_sort(&WP, EOM_TMP, qprs, 0, 7, "WabejDS_M");
    dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 0, 7, 2, 7, 0, SIJAB_lbl);
    dpd_buf4_axpy(&WP, &SIJAB, 1.0);
    dpd_buf4_close(&WP);
    dpd_buf4_init(&WM, EOM_TMP, C_irr, 0, 7, 0, 7, 0, "WabejDS_M");
    dpd_buf4_axpy(&WM, &SIJAB, -1.0);
    dpd_buf4_close(&WM);
    dpd_buf4_close(&SIJAB);

    /* Sijab += Wabej * Cie - Wabei * Cje */
    dpd_buf4_init(&WP, EOM_TMP, C_irr, 0, 7, 0, 7, 0, "WabejDS_P");
    dpd_buf4_init(&Wabei, CC_HBAR, H_IRR, 11, 7, 11, 7, 0, "Weiab");
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, Cme_lbl);
    dpd_contract244(&Cme, &Wabei, &WP, 1, 0, 0, 1.0, 0.0);
    dpd_file2_close(&Cme);
    dpd_buf4_close(&Wabei);
    dpd_buf4_sort(&WP, EOM_TMP, qprs, 0, 7, "WabejDS_M");
    dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 0, 7, 2, 7, 0, Sijab_lbl);
    dpd_buf4_axpy(&WP, &Sijab, 1.0);
    dpd_buf4_close(&WP);
    dpd_buf4_init(&WM, EOM_TMP, C_irr, 0, 7, 0, 7, 0, "WabejDS_M");
    dpd_buf4_axpy(&WM, &Sijab, -1.0);
    dpd_buf4_close(&WM);
    dpd_buf4_close(&Sijab);


    /* SIjAb += WAbEj * CIE - WAbeI * Cje */
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_buf4_init(&WAbEi, CC_HBAR, H_IRR, 11, 5, 11, 5, 0, "WEiAb");
    dpd_contract244(&CME, &WAbEi, &SIjAb, 1, 0, 0, 1.0, 1.0);
    dpd_file2_close(&CME);
    dpd_buf4_close(&WAbEi);
    dpd_buf4_init(&WaBeI, CC_HBAR, H_IRR, 10, 5, 10, 5, 0, "WeIaB (Ie,Ab)");
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, Cme_lbl);
    dpd_contract424(&WaBeI, &Cme, &SIjAb, 1, 1, 1, 1.0, 1.0);
    dpd_file2_close(&Cme);
    dpd_buf4_close(&WaBeI);
    dpd_buf4_close(&SIjAb);
  }

  else if (params.eom_ref == 2) {
    sprintf(CME_lbl, "%s %d", "CME", i);
    sprintf(Cme_lbl, "%s %d", "Cme", i);
    sprintf(SIJAB_lbl, "%s %d", "SIJAB", i);
    sprintf(Sijab_lbl, "%s %d", "Sijab", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    /* SIJAB += WABEJ * CIE - WABEI * CJE */
    dpd_buf4_init(&WP, EOM_TMP, C_irr, 0, 7, 0, 7, 0, "WabejDS_P");
    dpd_buf4_init(&WABEI, CC_HBAR, H_IRR, 21, 7, 21, 7, 0, "WEIAB");
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_contract244(&CME, &WABEI, &WP, 1, 0, 0, 1.0, 0.0);
    dpd_file2_close(&CME);
    dpd_buf4_close(&WABEI);
    dpd_buf4_sort(&WP, EOM_TMP, qprs, 0, 7, "WabejDS_M");
    dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 0, 7, 2, 7, 0, SIJAB_lbl);
    dpd_buf4_axpy(&WP, &SIJAB, 1.0);
    dpd_buf4_close(&WP);
    dpd_buf4_init(&WM, EOM_TMP, C_irr, 0, 7, 0, 7, 0, "WabejDS_M");
    dpd_buf4_axpy(&WM, &SIJAB, -1.0);
    dpd_buf4_close(&WM);
    dpd_buf4_close(&SIJAB);

    /* Sijab += Wabej * Cie - Wabei * Cje */
    dpd_buf4_init(&WP, EOM_TMP, C_irr, 10, 17, 10, 17, 0, "WabejDS_P");
    dpd_buf4_init(&Wabei, CC_HBAR, H_IRR, 31, 17, 31, 17, 0, "Weiab");
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 2, 3, Cme_lbl);
    dpd_contract244(&Cme, &Wabei, &WP, 1, 0, 0, 1.0, 0.0);
    dpd_file2_close(&Cme);
    dpd_buf4_close(&Wabei);
    dpd_buf4_sort(&WP, EOM_TMP, qprs, 10, 17, "WabejDS_M");
    dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 10, 17, 12, 17, 0, Sijab_lbl);
    dpd_buf4_axpy(&WP, &Sijab, 1.0);
    dpd_buf4_close(&WP);
    dpd_buf4_init(&WM, EOM_TMP, C_irr, 10, 17, 10, 17, 0, "WabejDS_M");
    dpd_buf4_axpy(&WM, &Sijab, -1.0);
    dpd_buf4_close(&WM);
    dpd_buf4_close(&Sijab);


    /* SIjAb += WAbEj * CIE - WAbeI * Cje */
    // start here
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 22, 28, 22, 28, 0, SIjAb_lbl);
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_buf4_init(&WAbEi, CC_HBAR, H_IRR, 26, 28, 26, 28, 0, "WEiAb");
    dpd_contract244(&CME, &WAbEi, &SIjAb, 1, 0, 0, 1.0, 1.0);
    dpd_file2_close(&CME);
    dpd_buf4_close(&WAbEi);
    dpd_buf4_init(&WaBeI, CC_HBAR, H_IRR, 24, 28, 24, 28, 0, "WeIaB (Ie,Ab)");
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 2, 3, Cme_lbl);
    dpd_contract424(&WaBeI, &Cme, &SIjAb, 1, 1, 1, 1.0, 1.0);
    dpd_file2_close(&Cme);
    dpd_buf4_close(&WaBeI);
    dpd_buf4_close(&SIjAb);
  }

#ifdef EOM_DEBUG
  check_sum("WabejDS",i,C_irr);
#endif
  return;
}

