#include <stdio.h>
#include <math.h>
#define EXTERN
#include "globals.h"

/* This function computes the H-bar doubles-doubles block contribution
from -0.5*P(ab)*Wnmfe*Cmnea*tijfb and +0.5*Wnmfe*Cimfe*tjnab to a
Sigma vector stored at Sigma plus 'i' */

void WmnefDD(int i, int irrep) {
  dpdbuf4 C2, T2, S2, D;
  dpdfile2 X;
  dpdbuf4 SIJAB, Sijab, SIjAb, B;
  dpdbuf4 CMNEF, Cmnef, CMnEf, F, tau;
  char CMNEF_lbl[32], Cmnef_lbl[32], CMnEf_lbl[32];
  char SIJAB_lbl[32], Sijab_lbl[32], SIjAb_lbl[32];
  int l,I,a,f,h,nirreps,*occpi,*virtpi,*openpi;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi; openpi = moinfo.openpi;

  sprintf(CMNEF_lbl, "%s %d", "CMNEF", i);
  sprintf(Cmnef_lbl, "%s %d", "Cmnef", i);
  sprintf(CMnEf_lbl, "%s %d", "CMnEf", i);
  sprintf(SIJAB_lbl, "%s %d", "SIJAB", i);
  sprintf(Sijab_lbl, "%s %d", "Sijab", i);
  sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

  /* XAF = CMNAE * WMNFE + CMnAe * WMnFe */
  /* SIJAB -= P(ab) XAF * TIJFB */
  dpd_file2_init(&X, EOM_TMP, irrep, 1, 1, "XFA");
  dpd_file2_mat_init(&X);
  dpd_file2_mat_wrt(&X);
  dpd_file2_mat_close(&X);
  dpd_buf4_init(&D, CC_DINTS, irrep, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");
  dpd_buf4_init(&CMNEF, EOM_CMNEF, irrep, 2, 5, 2, 7, 0, CMNEF_lbl);
  dpd_contract442(&D, &CMNEF, &X, 2, 2, 1.0, 1.0);
  dpd_buf4_close(&CMNEF);
  dpd_buf4_close(&D);
  dpd_buf4_init(&D, CC_DINTS, irrep, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_init(&CMnEf, EOM_CMnEf, irrep, 0, 5, 0, 5, 0, CMnEf_lbl);
  dpd_contract442(&D, &CMnEf, &X, 2, 2, 1.0, 1.0);
  dpd_buf4_close(&CMnEf);
  dpd_buf4_close(&D);

  dpd_buf4_init(&S2, EOM_TMP, irrep, 2, 5, 2, 5, 0, "SIJAB");
  dpd_buf4_init(&T2, CC_TAMPS, irrep, 2, 5, 2, 7, 0, "tIJAB");
  dpd_contract244(&X, &T2, &S2, 0, 2, 1, 1.0, 0.0);
  dpd_buf4_close(&T2);
  dpd_file2_close(&X);
  dpd_buf4_sort(&S2, EOM_TMP, pqsr, 2, 5, "SIJBA");
  dpd_buf4_init(&SIJAB, EOM_SIJAB, irrep, 2, 5, 2, 7, 0, SIJAB_lbl);
  dpd_buf4_axpy(&S2, &SIJAB, -1.0);
  dpd_buf4_close(&S2);
  dpd_buf4_init(&S2, EOM_TMP, irrep, 2, 5, 2, 5, 0, "SIJBA");
  dpd_buf4_axpy(&S2, &SIJAB, 1.0);
  dpd_buf4_close(&S2);
  dpd_buf4_close(&SIJAB);

  /* Xaf = Cmnae * Wmnfe + CmNaE * WmNfE */
  /* Sijab -= P(ab) Xfa * Tijfb */
  dpd_file2_init(&X, EOM_TMP, irrep, 1, 1, "Xfa");
  dpd_file2_mat_init(&X);
  dpd_file2_mat_wrt(&X);
  dpd_file2_mat_close(&X);
  dpd_buf4_init(&D, CC_DINTS, irrep, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");
  dpd_buf4_init(&Cmnef, EOM_Cmnef, irrep, 2, 5, 2, 7, 0, Cmnef_lbl);
  dpd_contract442(&D, &Cmnef, &X, 2, 2, 1.0, 1.0);
  dpd_buf4_close(&Cmnef);
  dpd_buf4_close(&D);

  dpd_buf4_init(&D, CC_DINTS, irrep, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_init(&CMnEf, EOM_TMP, irrep, 0, 5, 0, 5, 0, "CmNeF");
  dpd_contract442(&D, &CMnEf, &X, 2, 2, 1.0, 1.0);
  dpd_buf4_close(&CMnEf);
  dpd_buf4_close(&D);

  dpd_buf4_init(&S2, EOM_TMP, irrep, 2, 5, 2, 5, 0, "Sijab");
  dpd_buf4_init(&T2, CC_TAMPS, irrep, 2, 5, 2, 7, 0, "tijab");
  dpd_contract244(&X, &T2, &S2, 0, 2, 1, 1.0, 0.0);
  dpd_file2_close(&X);
  dpd_buf4_close(&T2);
  dpd_buf4_sort(&S2, EOM_TMP, pqsr, 2, 5, "Sijba");
  dpd_buf4_init(&Sijab, EOM_Sijab, irrep, 2, 5, 2, 7, 0, Sijab_lbl);
  dpd_buf4_axpy(&S2, &Sijab, -1.0);
  dpd_buf4_close(&S2);
  dpd_buf4_init(&S2, EOM_TMP, irrep, 2, 5, 2, 5, 0, "Sijba");
  dpd_buf4_axpy(&S2, &Sijab, 1.0);
  dpd_buf4_close(&S2);
  dpd_buf4_close(&Sijab);

  /* SIjAb += -XFA * TIjFb - TIjAf * Xfb */
  dpd_buf4_init(&SIjAb, EOM_SIjAb, irrep, 0, 5, 0, 5, 0, SIjAb_lbl);
  dpd_file2_init(&X, EOM_TMP, irrep, 1, 1, "XFA");
  dpd_buf4_init(&T2, CC_TAMPS, irrep, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract244(&X, &T2, &SIjAb, 0, 2, 1, -1.0, 1.0);
  dpd_buf4_close(&T2);
  dpd_file2_close(&X);
  dpd_file2_init(&X, EOM_TMP, irrep, 1, 1, "Xfa");
  dpd_buf4_init(&T2, CC_TAMPS, irrep, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract424(&T2, &X, &SIjAb, 3, 0, 0, -1.0, 1.0);
  dpd_buf4_close(&T2);
  dpd_file2_close(&X);
  dpd_buf4_close(&SIjAb);

#ifdef EOM_DEBUG
  check_sum("WmnefDD XAF",i,irrep);
#endif

  /* XLI = WLMEF * CIMEF + WLmEf * CImEf */
  /* SIJAB += P(IJ) XLI * TLJAB */
  dpd_file2_init(&X, EOM_TMP, irrep, 0, 0, "XLI");
  dpd_file2_mat_init(&X);
  dpd_file2_mat_wrt(&X);
  dpd_file2_mat_close(&X);
  dpd_buf4_init(&CMNEF, EOM_CMNEF, irrep, 0, 7, 2, 7, 0, CMNEF_lbl);
  dpd_buf4_init(&D, CC_DINTS, irrep, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
  dpd_contract442(&D, &CMNEF, &X, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&D);
  dpd_buf4_close(&CMNEF);
  dpd_buf4_init(&CMnEf, EOM_CMnEf, irrep, 0, 5, 0, 5, 0, CMnEf_lbl);
  dpd_buf4_init(&D, CC_DINTS, irrep, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_contract442(&D, &CMnEf, &X, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&D);
  dpd_buf4_close(&CMnEf);

  /*
    dpd_file2_mat_init(&X);
    dpd_file2_mat_rd(&X);
    for(h=0; h < nirreps; h++) {
    for(l=0; l<occpi[h]; l++)
    for(I=(occpi[h]-openpi[h]); I<occpi[h]; I++)
    X.matrix[h][l][I] = 0.0;
    for(I=0; I<occpi[h]; I++)
    for(l=(occpi[h]-openpi[h]); l<occpi[h]; l++)
    X.matrix[h][l][I] = 0.0;
    }
    dpd_file2_mat_wrt(&X);
    dpd_file2_mat_close(&X);
  */

  dpd_buf4_init(&S2, EOM_TMP, irrep, 0, 7, 0, 7, 0, "SIJAB");
  dpd_buf4_init(&T2, CC_TAMPS, irrep, 0, 7, 2, 7, 0, "tIJAB");
  dpd_contract244(&X, &T2, &S2, 0, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&T2);
  dpd_file2_close(&X);
  dpd_buf4_sort(&S2, EOM_TMP, qprs, 0, 7, "SJIAB");
  dpd_buf4_init(&SIJAB, EOM_SIJAB, irrep, 0, 7, 2, 7, 0, SIJAB_lbl);
  dpd_buf4_axpy(&S2, &SIJAB, -1.0);
  dpd_buf4_close(&S2);
  dpd_buf4_init(&S2, EOM_TMP, irrep, 0, 7, 0, 7, 0, "SJIAB");
  dpd_buf4_axpy(&S2, &SIJAB, 1.0);
  dpd_buf4_close(&S2);
  dpd_buf4_close(&SIJAB);


  /* Xli = Wlmef * Cimef + WlMeF * CiMeF */
  /* Sijab += P(ij) Xli * Tljab */
  dpd_file2_init(&X, EOM_TMP, irrep, 0, 0, "Xli");
  dpd_file2_mat_init(&X);
  dpd_file2_mat_wrt(&X);
  dpd_file2_mat_close(&X);
  dpd_buf4_init(&Cmnef, EOM_Cmnef, irrep, 0, 7, 2, 7, 0, Cmnef_lbl);
  dpd_buf4_init(&D, CC_DINTS, irrep, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
  dpd_contract442(&D, &Cmnef, &X, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&D);
  dpd_buf4_close(&Cmnef);
  dpd_buf4_init(&CMnEf, EOM_TMP, irrep, 0, 5, 0, 5, 0, "CmNeF");
  dpd_buf4_init(&D, CC_DINTS, irrep, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_contract442(&D, &CMnEf, &X, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&D);
  dpd_buf4_close(&CMnEf);
  dpd_buf4_init(&S2, EOM_TMP, irrep, 0, 7, 0, 7, 0, "Sijab");
  dpd_buf4_init(&T2, CC_TAMPS, irrep, 0, 7, 2, 7, 0, "tijab");
  dpd_contract244(&X, &T2, &S2, 0, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&T2);
  dpd_file2_close(&X);
  dpd_buf4_sort(&S2, EOM_TMP, qprs, 0, 7, "Sjiab");
  dpd_buf4_init(&Sijab, EOM_Sijab, irrep, 0, 7, 2, 7, 0, Sijab_lbl);
  dpd_buf4_axpy(&S2, &Sijab, -1.0);
  dpd_buf4_close(&S2);
  dpd_buf4_init(&S2, EOM_TMP, irrep, 0, 7, 0, 7, 0, "Sjiab");
  dpd_buf4_axpy(&S2, &Sijab, 1.0);
  dpd_buf4_close(&S2);
  dpd_buf4_close(&Sijab);

  /* SIjAb += -XLI * TLjAb - Xli * TIlAb */
  dpd_buf4_init(&SIjAb, EOM_SIjAb, irrep, 0, 5, 0, 5, 0, SIjAb_lbl);
  dpd_file2_init(&X, EOM_TMP, irrep, 0, 0, "XLI");
  dpd_buf4_init(&T2, CC_TAMPS, irrep, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract244(&X, &T2, &SIjAb, 0, 0, 0, -1.0, 1.0);
  dpd_file2_close(&X);
  dpd_file2_init(&X, EOM_TMP, irrep, 0, 0, "Xli");
  dpd_contract424(&T2, &X, &SIjAb, 1, 0, 1, -1.0, 1.0);
  dpd_buf4_close(&T2);
  dpd_file2_close(&X);
  dpd_buf4_close(&SIjAb);

#ifdef EOM_DEBUG
  check_sum("WmnefDD XLI",i,irrep);
#endif

  return;
}
