#include <stdio.h>
#include <math.h>
#define EXTERN
#include "globals.h"

/* This function computes the H-bar doubles-singles block contribution
of Wnmje to a Sigma vector stored at Sigma plus 'i' */

void WnmjeDS(int i, int irrep) {
  dpdfile2 CME, Cme, XNJ, Xnj;
  dpdbuf4 SIJAB, Sijab, SIjAb;
  dpdbuf4 WMNIE, Wmnie, WMnIe, WmNiE, WM, WP;
  dpdbuf4 TIJAB, TIjAb, Tijab;
  char CME_lbl[32], Cme_lbl[32], SIJAB_lbl[32], Sijab_lbl[32], SIjAb_lbl[32];

  sprintf(CME_lbl, "%s %d", "CME", i);
  sprintf(Cme_lbl, "%s %d", "Cme", i);
  sprintf(SIJAB_lbl, "%s %d", "SIJAB", i);
  sprintf(Sijab_lbl, "%s %d", "Sijab", i);
  sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

  /* Form XNJ intermediates */
  /* XNJ = CME * WNMJE + Cme * WNmJe */
  dpd_file2_init(&XNJ, EOM_TMP, irrep, 0, 0, "XNJ");
/*
  dpd_file2_mat_init(&XNJ);
  dpd_file2_mat_wrt(&XNJ);
*/
  dpd_file2_scm(&XNJ, 0.0);
  dpd_file2_init(&CME, EOM_CME, irrep, 0, 1, CME_lbl);
  dpd_buf4_init(&WMNIE, CC_HBAR, irrep, 0, 11, 2, 11, 0, "WMNIE");
  dpd_dot23(&CME, &WMNIE, &XNJ, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&WMNIE);
  dpd_file2_close(&CME);
  dpd_file2_init(&Cme, EOM_Cme, irrep, 0, 1, Cme_lbl);
  dpd_buf4_init(&WMnIe, CC_HBAR, irrep, 0, 11, 0, 11, 0, "WMnIe");
  dpd_dot23(&Cme, &WMnIe, &XNJ, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&WMnIe);
  dpd_file2_close(&Cme);
  dpd_file2_close(&XNJ);


  /* Xnj = Cme * Wnmje + CME * WnMjE */
  dpd_file2_init(&Xnj, EOM_TMP, irrep, 0, 0, "Xnj");
/*
  dpd_file2_mat_init(&Xnj);
  dpd_file2_mat_wrt(&Xnj);
*/
  dpd_file2_scm(&Xnj, 0.0);
  dpd_file2_init(&Cme, EOM_Cme, irrep, 0, 1, Cme_lbl);
  dpd_buf4_init(&Wmnie, CC_HBAR, irrep, 0, 11, 2, 11, 0, "Wmnie");
  dpd_dot23(&Cme, &Wmnie, &Xnj, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&Wmnie);
  dpd_file2_close(&Cme);
  dpd_file2_init(&CME, EOM_CME, irrep, 0, 1, CME_lbl);
  dpd_buf4_init(&WmNiE, CC_HBAR, irrep, 0, 11, 0, 11, 0, "WmNiE");
  dpd_dot23(&CME, &WmNiE, &Xnj, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&WmNiE);
  dpd_file2_close(&CME);
  dpd_file2_close(&Xnj);

  /* SIJAB -= XNJ * TINAB + XNI * TJNAB */
  dpd_buf4_init(&WM, EOM_TMP, irrep, 0, 7, 0, 7, 0, "WnmjeDS_M");
  dpd_buf4_init(&TIJAB, CC_TAMPS, irrep, 0, 7, 2, 7, 0, "tIJAB");
  dpd_file2_init(&XNJ, EOM_TMP, irrep, 0, 0, "XNJ");
  dpd_contract424(&TIJAB, &XNJ, &WM, 1, 0, 1, 1.0, 0.0);
  dpd_file2_close(&XNJ);
  dpd_buf4_close(&TIJAB);
  dpd_buf4_sort(&WM, EOM_TMP, qprs, 0, 7, "WnmjeDS_P");
  dpd_buf4_init(&SIJAB, EOM_SIJAB, irrep, 0, 7, 2, 7, 0, SIJAB_lbl);
  dpd_buf4_axpy(&WM, &SIJAB, -1.0);
  dpd_buf4_close(&WM);
  dpd_buf4_init(&WP, EOM_TMP, irrep, 0, 7, 0, 7, 0, "WnmjeDS_P");
  dpd_buf4_axpy(&WP, &SIJAB, 1.0);
  dpd_buf4_close(&WP);
  dpd_buf4_close(&SIJAB);

  /* Sijab -= Xnj * Tinab + Xni * Tjnab */
  dpd_buf4_init(&WM, EOM_TMP, irrep, 0, 7, 0, 7, 0, "WnmjeDS_M");
  dpd_buf4_init(&Tijab, CC_TAMPS, irrep, 0, 7, 2, 7, 0, "tijab");
  dpd_file2_init(&Xnj, EOM_TMP, irrep, 0, 0, "Xnj");
  dpd_contract424(&Tijab, &Xnj, &WM, 1, 0, 1, 1.0, 0.0);
  dpd_file2_close(&Xnj);
  dpd_buf4_close(&Tijab);
  dpd_buf4_sort(&WM, EOM_TMP, qprs, 0, 7, "WnmjeDS_P");
  dpd_buf4_init(&Sijab, EOM_Sijab, irrep, 0, 7, 2, 7, 0, Sijab_lbl);
  dpd_buf4_axpy(&WM, &Sijab, -1.0);
  dpd_buf4_close(&WM);
  dpd_buf4_init(&WP, EOM_TMP, irrep, 0, 7, 0, 7, 0, "WnmjeDS_P");
  dpd_buf4_axpy(&WP, &Sijab, 1.0);
  dpd_buf4_close(&WP);
  dpd_buf4_close(&Sijab);

  /* SIjAb -= Xnj * tInAb + XNI * TjNAb */
  dpd_buf4_init(&SIjAb, EOM_SIjAb, irrep, 0, 5, 0, 5, 0, SIjAb_lbl);
  dpd_buf4_init(&TIjAb, CC_TAMPS, irrep, 0, 5, 0, 5, 0, "tIjAb");
  dpd_file2_init(&Xnj, EOM_TMP, irrep, 0, 0, "Xnj");
  dpd_contract424(&TIjAb, &Xnj, &SIjAb, 1, 0, 1, -1.0, 1.0);
  dpd_file2_close(&Xnj);
  dpd_file2_init(&XNJ, EOM_TMP, irrep, 0, 0, "XNJ");
  dpd_contract244(&XNJ, &TIjAb, &SIjAb, 0, 0, 0, -1.0, 1.0);
  dpd_file2_close(&XNJ);
  dpd_buf4_close(&TIjAb);
  dpd_buf4_close(&SIjAb);

#ifdef EOM_DEBUG
  check_sum("WnmjeDS",i,irrep);
#endif

  return;
}
