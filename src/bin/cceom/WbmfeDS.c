#include <stdio.h>
#include <math.h>
#define EXTERN
#include "globals.h"

/* This function computes the H-bar doubles-singles block contribution
   of Wbmfe to a Sigma vector stored at Sigma plus 'i' */

void WbmfeDS(int i, int C_irr) {
  dpdfile2 CME, Cme, XBF, Xbf;
  dpdbuf4 SIJAB, Sijab, SIjAb;
  dpdbuf4 WAMEF, Wamef, WAmEf, WaMeF, WM, WP, W, Z;
  dpdbuf4 TIJAB, TIjAb, Tijab;
  char CME_lbl[32], Cme_lbl[32], SIJAB_lbl[32], Sijab_lbl[32], SIjAb_lbl[32];

  if (params.eom_ref == 0) { // RHF
    sprintf(CME_lbl, "%s %d", "CME", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    /* Form Xbf intermediates */
    dpd_file2_init(&XBF, EOM_TMP, C_irr, 1, 1, "XBF");
    dpd_file2_scm(&XBF, 0.0);
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 10, 5, 10, 5, 0, "WAmEf 2(mA,Ef) - (mA,fE)");
    dpd_dot14(&CME, &W, &XBF, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_file2_close(&CME);

    dpd_buf4_init(&Z, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "WbmfeDS Z(Ij,Ab)");
    dpd_buf4_init(&TIjAb, CC_TAMPS, H_IRR, 0, 5, 0, 5, 0, "tIjAb");
    dpd_contract424(&TIjAb, &XBF, &Z, 3, 1, 0, 1.0, 0.0);
    dpd_buf4_close(&TIjAb);
    dpd_buf4_sort(&Z, EOM_TMP, qpsr, 0, 5, "WbmfeDS Z(jI,bA)");

    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    dpd_buf4_axpy(&Z, &SIjAb, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "WbmfeDS Z(jI,bA)");
    dpd_buf4_axpy(&Z, &SIjAb, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&SIjAb);
    dpd_file2_close(&XBF);
  }

  else { // ROHF
    sprintf(CME_lbl, "%s %d", "CME", i);
    sprintf(Cme_lbl, "%s %d", "Cme", i);
    sprintf(SIJAB_lbl, "%s %d", "SIJAB", i);
    sprintf(Sijab_lbl, "%s %d", "Sijab", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    /* Form Xbf intermediates */
    /* XBF = CME * WBMFE + Cme * WBmFe */
    dpd_file2_init(&XBF, EOM_TMP, C_irr, 1, 1, "XBF");
    dpd_file2_scm(&XBF, 0.0);
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_buf4_init(&WAMEF, CC_HBAR, H_IRR, 10, 5, 10, 7, 0, "WAMEF");
    dpd_dot14(&CME, &WAMEF, &XBF, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&WAMEF);
    dpd_file2_close(&CME);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, Cme_lbl);
    dpd_buf4_init(&WAmEf, CC_HBAR, H_IRR, 10, 5, 10, 5, 0, "WAmEf");
    dpd_dot14(&Cme, &WAmEf, &XBF, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&WAmEf);
    dpd_file2_close(&Cme);
    dpd_file2_close(&XBF);

    /* Xbf = Cme * Wbmfe + CME * WbMfE */
    dpd_file2_init(&Xbf, EOM_TMP, C_irr, 1, 1, "Xbf");
    dpd_file2_scm(&Xbf, 0.0);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, Cme_lbl);
    dpd_buf4_init(&Wamef, CC_HBAR, H_IRR, 10, 5, 10, 7, 0, "Wamef");
    dpd_dot14(&Cme, &Wamef, &Xbf, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&Wamef);
    dpd_file2_close(&Cme);
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_buf4_init(&WaMeF, CC_HBAR, H_IRR, 10, 5, 10, 5, 0, "WaMeF");
    dpd_dot14(&CME, &WaMeF, &Xbf, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&WaMeF);
    dpd_file2_close(&CME);
    dpd_file2_close(&Xbf);

    /* SIJAB += XBF * TIJAF - XAF * TIJBF */
    dpd_buf4_init(&WP, EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WbmfeDS_P");
    dpd_file2_init(&XBF, EOM_TMP, C_irr, 1, 1, "XBF");
    dpd_buf4_init(&TIJAB, CC_TAMPS, H_IRR, 2, 5, 2, 7, 0, "tIJAB");
    dpd_contract424(&TIJAB, &XBF, &WP, 3, 1, 0, 1.0, 0.0);
    dpd_buf4_close(&TIJAB);
    dpd_file2_close(&XBF);
    dpd_buf4_sort(&WP, EOM_TMP, pqsr, 2, 5, "WbmfeDS_M"); 
    dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 2, 5, 2, 7, 0, SIJAB_lbl);
    dpd_buf4_axpy(&WP, &SIJAB, 1.0);
    dpd_buf4_close(&WP);
    dpd_buf4_init(&WM, EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WbmfeDS_M");
    dpd_buf4_axpy(&WM, &SIJAB, -1.0);
    dpd_buf4_close(&WM);
    dpd_buf4_close(&SIJAB);

    /* Sijab += Xbf * Tijaf - Xaf * Tijbf */
    dpd_buf4_init(&WP, EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WbmfeDS_P");
    dpd_file2_init(&Xbf, EOM_TMP, C_irr, 1, 1, "Xbf");
    dpd_buf4_init(&Tijab, CC_TAMPS, H_IRR, 2, 5, 2, 7, 0, "tijab");
    dpd_contract424(&Tijab, &Xbf, &WP, 3, 1, 0, 1.0, 0.0);
    dpd_buf4_close(&Tijab);
    dpd_file2_close(&Xbf);
    dpd_buf4_sort(&WP, EOM_TMP, pqsr, 2, 5, "WbmfeDS_M");
    dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 2, 5, 2, 7, 0, Sijab_lbl);
    dpd_buf4_axpy(&WP, &Sijab, 1.0);
    dpd_buf4_close(&WP);
    dpd_buf4_init(&WM, EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WbmfeDS_M");
    dpd_buf4_axpy(&WM, &Sijab, -1.0);
    dpd_buf4_close(&WM);
    dpd_buf4_close(&Sijab);

    /* SIjAb += Xbf * tIjAf + XAF * TIjbF */
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    dpd_buf4_init(&TIjAb, CC_TAMPS, H_IRR, 0, 5, 0, 5, 0, "tIjAb");
    dpd_file2_init(&Xbf, EOM_TMP, C_irr, 1, 1, "Xbf");
    dpd_contract424(&TIjAb, &Xbf, &SIjAb, 3, 1, 0, 1.0, 1.0);
    dpd_file2_close(&Xbf);
    dpd_file2_init(&XBF, EOM_TMP, C_irr, 1, 1, "XBF");
    dpd_contract244(&XBF, &TIjAb, &SIjAb, 1, 2, 1, 1.0, 1.0);
    dpd_file2_close(&XBF);
    dpd_buf4_close(&TIjAb);
    dpd_buf4_close(&SIjAb);
  }

#ifdef EOM_DEBUG
  check_sum("WbmfeDS",i,C_irr);
#endif
  return;
}
