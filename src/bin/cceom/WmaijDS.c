#include <stdio.h>
#include <math.h>
#define EXTERN
#include "globals.h"

/* This function computes the H-bar doubles-singles block contribution
of Wmaij to a Sigma vector stored at Sigma plus 'i' */

void WmaijDS(int i, int C_irr) {
  dpdfile2 CME, Cme;
  dpdbuf4 SIJAB, Sijab, SIjAb;
  dpdbuf4 WP, WM, WMBIJ, Wmbij, WMbIj, WmBiJ;
  char CME_lbl[32], Cme_lbl[32], SIJAB_lbl[32], Sijab_lbl[32], SIjAb_lbl[32];

  sprintf(CME_lbl, "%s %d", "CME", i);
  sprintf(Cme_lbl, "%s %d", "Cme", i);
  sprintf(SIJAB_lbl, "%s %d", "SIJAB", i);
  sprintf(Sijab_lbl, "%s %d", "Sijab", i);
  sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

 /* SIJAB += WMAIJ * CMB - WMBIJ * CMA */
  dpd_buf4_init(&WM, EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WmaijDS_M");
  dpd_buf4_init(&WMBIJ, CC_HBAR, H_IRR, 10, 2, 10, 2, 0, "WMBIJ");
  dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
  dpd_contract244(&CME, &WMBIJ, &WM, 0, 0, 1, 1.0, 0.0);
  dpd_file2_close(&CME);
  dpd_buf4_close(&WMBIJ);
  dpd_buf4_sort(&WM, EOM_TMP, pqsr, 2, 5, "WmaijDS_P");
  dpd_buf4_init(&WP, EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WmaijDS_P");
  dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 2, 5, 2, 7, 0, SIJAB_lbl);
  dpd_buf4_axpy(&WM, &SIJAB, -1.0);
  dpd_buf4_axpy(&WP, &SIJAB,  1.0);
  dpd_buf4_close(&SIJAB);
  dpd_buf4_close(&WM);
  dpd_buf4_close(&WP);


 /* Sijab += Wmaij * Cmb - Wmbij * Cma */
  dpd_buf4_init(&WM, EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WmaijDS_M");
  dpd_buf4_init(&Wmbij, CC_HBAR, H_IRR, 10, 2, 10, 2, 0, "Wmbij");
  dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, Cme_lbl);
  dpd_contract244(&Cme, &Wmbij, &WM, 0, 0, 1, 1.0, 0.0);
  dpd_file2_close(&Cme);
  dpd_buf4_close(&Wmbij);
  dpd_buf4_sort(&WM, EOM_TMP, pqsr, 2, 5, "WmaijDS_P");
  dpd_buf4_init(&WP, EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WmaijDS_P");
  dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 2, 5, 2, 7, 0, Sijab_lbl);
  dpd_buf4_axpy(&WM, &Sijab, -1.0);
  dpd_buf4_axpy(&WP, &Sijab,  1.0);
  dpd_buf4_close(&Sijab);
  dpd_buf4_close(&WM);
  dpd_buf4_close(&WP);


 /* SIjAb += WmAIj * Cmb - WMbIj * CMA */
  dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
  dpd_buf4_init(&WmBiJ, CC_HBAR, H_IRR, 11, 0, 11, 0, 0, "WmBiJ (Bm,Ji)");
  dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, Cme_lbl);
  dpd_contract424(&WmBiJ, &Cme, &SIjAb, 1, 0, 0, -1.0, 0.0);
  dpd_file2_close(&Cme);
  dpd_buf4_close(&WmBiJ);

  dpd_buf4_init(&WMbIj, CC_HBAR, H_IRR, 10, 0, 10, 0, 0, "WMbIj");
  dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
  dpd_contract244(&CME, &WMbIj, &SIjAb, 0, 0, 1, -1.0, 1.0);
  dpd_file2_close(&CME);
  dpd_buf4_close(&WMbIj);
  dpd_buf4_close(&SIjAb);

#ifdef EOM_DEBUG
  check_sum("WmaijDS",i,C_irr);
#endif

  return;
}
