#include <stdio.h>
#include <math.h>
#define EXTERN
#include "globals.h"

/* This function computes the H-bar doubles-doubles block contribution
to a Sigma vector stored at Sigma plus 'i' */

void FDD(int i, int irrep) {
  dpdfile2 FAE, Fae, FMI, Fmi;
  dpdbuf4 SIJAB, Sijab, SIjAb, FP, FM;
  dpdbuf4 CMNEF, Cmnef, CMnEf;
  char CMNEF_lbl[32], Cmnef_lbl[32], CMnEf_lbl[32], CmNeF_lbl[32];
  char SIJAB_lbl[32], Sijab_lbl[32], SIjAb_lbl[32];

  sprintf(CMNEF_lbl, "%s %d", "CMNEF", i);
  sprintf(Cmnef_lbl, "%s %d", "Cmnef", i);
  sprintf(CMnEf_lbl, "%s %d", "CMnEf", i);
  sprintf(CmNeF_lbl, "%s %d", "CmNeF", i);
  sprintf(SIJAB_lbl, "%s %d", "SIJAB", i);
  sprintf(Sijab_lbl, "%s %d", "Sijab", i);
  sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

  /* SIJAB += FBE * CIJAE - FAE * CIJBE */
  dpd_buf4_init(&FP, EOM_TMP, irrep, 2, 5, 2, 5, 0, "FDD_FBEP");
  dpd_buf4_init(&CMNEF, EOM_CMNEF, irrep, 2, 5, 2, 7, 0, CMNEF_lbl);
  dpd_file2_init(&FAE, CC_OEI, irrep, 1, 1, "FAE");
  dpd_contract424(&CMNEF, &FAE, &FP, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&FAE);
  dpd_buf4_close(&CMNEF);
  dpd_buf4_sort(&FP, EOM_TMP, pqsr, 2, 5, "FDD_FBEM");
  dpd_buf4_init(&SIJAB, EOM_SIJAB, irrep, 2, 5, 2, 7, 0, SIJAB_lbl);
  dpd_buf4_axpy(&FP, &SIJAB, 1.0);
  dpd_buf4_close(&FP);
  dpd_buf4_init(&FM, EOM_TMP, irrep, 2, 5, 2, 5, 0, "FDD_FBEM");
  dpd_buf4_axpy(&FM, &SIJAB, -1.0);
  dpd_buf4_close(&FM);
  dpd_buf4_close(&SIJAB);

  /* Sijab += Fbe * Cijae - Fae * Cijbe */
  dpd_buf4_init(&FP, EOM_TMP, irrep, 2, 5, 2, 5, 0, "FDD_FBEP");
  dpd_buf4_init(&Cmnef, EOM_Cmnef, irrep, 2, 5, 2, 7, 0, Cmnef_lbl);
  dpd_file2_init(&Fae, CC_OEI, irrep, 1, 1, "Fae");
  dpd_contract424(&Cmnef, &Fae, &FP, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&Fae);
  dpd_buf4_close(&Cmnef);
  dpd_buf4_sort(&FP, EOM_TMP, pqsr, 2, 5, "FDD_FBEM");
  dpd_buf4_init(&Sijab, EOM_Sijab, irrep, 2, 5, 2, 7, 0, Sijab_lbl);
  dpd_buf4_axpy(&FP, &Sijab, 1.0);
  dpd_buf4_close(&FP);
  dpd_buf4_init(&FM, EOM_TMP, irrep, 2, 5, 2, 5, 0, "FDD_FBEM");
  dpd_buf4_axpy(&FM, &Sijab, -1.0);
  dpd_buf4_close(&FM);
  dpd_buf4_close(&Sijab);

  /* SIjAb += Fbe * CIjAe - FAE * CIjbE */
  dpd_buf4_init(&SIjAb, EOM_SIjAb, irrep, 0, 5, 0, 5, 0, SIjAb_lbl);
  dpd_buf4_init(&CMnEf, EOM_CMnEf, irrep, 0, 5, 0, 5, 0, CMnEf_lbl);
  dpd_file2_init(&Fae, CC_OEI, irrep, 1, 1, "Fae");
  dpd_contract424(&CMnEf, &Fae, &SIjAb, 3, 1, 0, 1.0, 1.0);
  dpd_file2_close(&Fae);
  dpd_file2_init(&FAE, CC_OEI, irrep, 1, 1, "FAE");
  dpd_contract244(&FAE, &CMnEf, &SIjAb, 1, 2, 1, 1.0, 1.0);
  dpd_file2_close(&FAE);
  dpd_buf4_close(&CMnEf);
  dpd_buf4_close(&SIjAb);

#ifdef EOM_DEBUG
  check_sum("FDD_Fbe",i,irrep);
#endif

  /* SIJAB -= FMJ * CIMAB - FMI * CJMAB */
  dpd_buf4_init(&FM, EOM_TMP, irrep, 0, 7, 0, 7, 0, "FDD_FMJM");
  dpd_buf4_init(&CMNEF, EOM_CMNEF, irrep, 0, 7, 2, 7, 0, CMNEF_lbl);
  dpd_file2_init(&FMI, CC_OEI, irrep, 0, 0, "FMI");
  dpd_contract424(&CMNEF, &FMI, &FM, 1, 0, 1, 1.0, 0.0);
  dpd_file2_close(&FMI);
  dpd_buf4_close(&CMNEF);
  dpd_buf4_sort(&FM, EOM_TMP, qprs, 0, 7, "FDD_FMJP");
  dpd_buf4_init(&SIJAB, EOM_SIJAB, irrep, 0, 7, 2, 7, 0, SIJAB_lbl);
  dpd_buf4_axpy(&FM, &SIJAB, -1.0);
  dpd_buf4_close(&FM);
  dpd_buf4_init(&FP, EOM_TMP, irrep, 0, 7, 0, 7, 0, "FDD_FMJP");
  dpd_buf4_axpy(&FP, &SIJAB, 1.0);
  dpd_buf4_close(&FP);
  dpd_buf4_close(&SIJAB);

  /* Sijab -= Fmj * Cimab - Fmi * Cjmab */
  dpd_buf4_init(&FM, EOM_TMP, irrep, 0, 7, 0, 7, 0, "FDD_FMJM");
  dpd_buf4_init(&Cmnef, EOM_Cmnef, irrep, 0, 7, 2, 7, 0, Cmnef_lbl);
  dpd_file2_init(&Fmi, CC_OEI, irrep, 0, 0, "Fmi");
  dpd_contract424(&Cmnef, &Fmi, &FM, 1, 0, 1, 1.0, 0.0);
  dpd_file2_close(&Fmi);
  dpd_buf4_close(&Cmnef);
  dpd_buf4_sort(&FM, EOM_TMP, qprs, 0, 7, "FDD_FMJP");
  dpd_buf4_init(&Sijab, EOM_Sijab, irrep, 0, 7, 2, 7, 0, Sijab_lbl);
  dpd_buf4_axpy(&FM, &Sijab, -1.0);
  dpd_buf4_close(&FM);
  dpd_buf4_init(&FP, EOM_TMP, irrep, 0, 7, 0, 7, 0, "FDD_FMJP");
  dpd_buf4_axpy(&FP, &Sijab, 1.0);
  dpd_buf4_close(&FP);
  dpd_buf4_close(&Sijab);

  /* SIjAb -= Fmj * CImAb - FMI * CjMAb */
  dpd_buf4_init(&SIjAb, EOM_SIjAb, irrep, 0, 5, 0, 5, 0, SIjAb_lbl);
  dpd_buf4_init(&CMnEf, EOM_CMnEf, irrep, 0, 5, 0, 5, 0, CMnEf_lbl);
  dpd_file2_init(&Fmi, CC_OEI, irrep, 0, 0, "Fmi");
  dpd_contract424(&CMnEf, &Fmi, &SIjAb, 1, 0, 1, -1.0, 1.0);
  dpd_file2_close(&Fmi);
  dpd_file2_init(&FMI, CC_OEI, irrep, 0, 0, "FMI");
  dpd_contract244(&FMI, &CMnEf, &SIjAb, 0, 0, 0, -1.0, 1.0);
  dpd_file2_close(&FMI);
  dpd_buf4_close(&CMnEf);
  dpd_buf4_close(&SIjAb);

#ifdef EOM_DEBUG
  check_sum("FDD_Fmj",i,irrep);
#endif

  return;
}
