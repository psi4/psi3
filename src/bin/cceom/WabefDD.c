#include <stdio.h>
#include <math.h>
#define EXTERN
#include "globals.h"

/* This function computes the H-bar doubles-doubles block contribution
   from Wabef to a Sigma vector stored at Sigma plus 'i' */

void WabefDD(int i, int C_irr) {
  dpdfile2 tIA, tia, SIA, Sia;
  dpdbuf4 SIJAB, Sijab, SIjAb, B;
  dpdbuf4 CMNEF, Cmnef, CMnEf, X, F, tau, D, WM, WP, Z;
  char CMNEF_lbl[32], Cmnef_lbl[32], CMnEf_lbl[32];
  char SIJAB_lbl[32], Sijab_lbl[32], SIjAb_lbl[32], SIA_lbl[32], Sia_lbl[32];

  if (params.eom_ref == 0) { // RHF
    /* SIjAb += WAbEf*CIjEf */
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);
    sprintf(CMnEf_lbl, "%s %d", "CMnEf", i);

    /* SIjAb += <Ab|Ef> CIjEf -- allow out of core algorithm */

#ifdef TIME_CCEOM
    timer_on("WabefDD Z");
#endif

    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, CMnEf_lbl);
    dpd_buf4_init(&Z, EOM_TMP, C_irr, 5, 0, 5, 0, 0, "WabefDD Z(Ab,Ij)");
    dpd_buf4_init(&B, CC_BINTS, H_IRR, 5, 5, 5, 5, 0, "B <ab|cd>");
    dpd_contract444(&B, &CMnEf, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&B);
    dpd_buf4_close(&CMnEf);
    dpd_buf4_sort(&Z, EOM_TMP, rspq, 0, 5, "WabefDD Z(Ij,Ab)");
    dpd_buf4_close(&Z);

    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    dpd_buf4_init(&Z, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "WabefDD Z(Ij,Ab)");
    dpd_buf4_axpy(&Z, &SIjAb, 1);
    dpd_buf4_close(&Z);

#ifdef TIME_CCEOM
    timer_off("WabefDD Z");
#endif

    /* construct XIjMb = CIjEf * <mb|ef> */
    dpd_buf4_init(&X, EOM_TMP, C_irr, 0, 10, 0, 10, 0, "WabefDD X(Ij,Mb)");
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, CMnEf_lbl);
    dpd_buf4_init(&F, CC_FINTS, H_IRR, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_contract444(&CMnEf, &F, &X, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&F);
    dpd_buf4_close(&CMnEf);

    dpd_buf4_init(&Z, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "WabefDD Z(Ij,Ab)");
    dpd_file2_init(&tIA, CC_OEI, H_IRR, 0, 1, "tIA");
//fprintf(outfile,"\n begin contract244 in WabefDD\n");
    dpd_contract244(&tIA, &X, &Z, 0, 2, 1, 1.0, 0.0);
    dpd_file2_close(&tIA);
//dpd_buf4_print(&Z,outfile,1);
    dpd_buf4_close(&X);

    dpd_buf4_sort(&Z, EOM_TMP, qpsr, 0, 5, "WabefDD Z(jI,bA)");
    dpd_buf4_axpy(&Z, &SIjAb, -1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "WabefDD Z(jI,bA)");
    dpd_buf4_axpy(&Z, &SIjAb, -1.0);
    dpd_buf4_close(&Z);

    /* SIjAb += tau_MnAb <Mn||ef> CIjEf */
    dpd_buf4_init(&X, EOM_TMP, C_irr, 0, 0, 0, 0, 0, "WabefDD XIjMn");
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, CMnEf_lbl);
    dpd_buf4_init(&D, CC_DINTS, H_IRR, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract444(&CMnEf, &D, &X, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&CMnEf);
    dpd_buf4_init(&tau, CC_TAMPS, H_IRR, 0, 5, 0, 5, 0, "tauIjAb");
    dpd_contract444(&X, &tau, &SIjAb, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&tau);
    dpd_buf4_close(&X);
    dpd_buf4_close(&SIjAb);
  }

  else { // ROHF
    sprintf(CMNEF_lbl, "%s %d", "CMNEF", i);
    sprintf(Cmnef_lbl, "%s %d", "Cmnef", i);
    sprintf(SIJAB_lbl, "%s %d", "SIJAB", i);
    sprintf(Sijab_lbl, "%s %d", "Sijab", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    /* SIJAB += WABEF*CIJEF */
    dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, CMNEF_lbl);
    dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 2, 7, 2, 7, 0, SIJAB_lbl);
    dpd_buf4_init(&B, CC_BINTS, H_IRR, 7, 7, 5, 5, 1, "B <ab|cd>");
    dpd_contract444(&CMNEF, &B, &SIJAB, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&B);
    dpd_buf4_close(&SIJAB);
    dpd_buf4_init(&X, EOM_TMP, C_irr, 2, 10, 2, 10, 0, "XIJMA");
    dpd_buf4_init(&F, CC_FINTS, H_IRR, 10, 7, 10, 5, 1, "F <ia|bc>");
    dpd_contract444(&CMNEF, &F, &X, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&F);
    dpd_buf4_close(&CMNEF);
    dpd_buf4_init(&WM, EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WabefDD_M");
    dpd_file2_init(&tIA, CC_OEI, H_IRR, 0, 1, "tIA");
    dpd_contract244(&tIA, &X, &WM, 0, 2, 1, 1.0, 0.0);
    dpd_file2_close(&tIA);
    dpd_buf4_sort(&WM, EOM_TMP, pqsr, 2, 5, "WabefDD_P");
    dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 2, 5, 2, 7, 0, SIJAB_lbl);
    dpd_buf4_axpy(&WM, &SIJAB, -1.0);
    dpd_buf4_close(&WM);
    dpd_buf4_init(&WP, EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WabefDD_P");
    dpd_buf4_axpy(&WP, &SIJAB, 1.0);
    dpd_buf4_close(&WP);
    dpd_buf4_close(&SIJAB);
    dpd_buf4_init(&X, EOM_TMP, C_irr, 2, 2, 2, 2, 0, "XIJMN");
    dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, CMNEF_lbl);
    dpd_buf4_init(&D, CC_DINTS, H_IRR, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
    dpd_contract444(&CMNEF, &D, &X, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&CMNEF);
    dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 2, 7, 2, 7, 0, SIJAB_lbl);
    dpd_buf4_init(&tau, CC_TAMPS, H_IRR, 2, 7, 2, 7, 0, "tauIJAB");
    dpd_contract444(&X, &tau, &SIJAB, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&tau);
    dpd_buf4_close(&SIJAB);
    dpd_buf4_close(&X);

#ifdef EOM_DEBUG
  check_sum("WabefDD",i,C_irr);
#endif

    /* Sijab += Wabef*Cijef */
    dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 2, 7, 2, 7, 0, Cmnef_lbl);
    dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 2, 7, 2, 7, 0, Sijab_lbl);
    dpd_buf4_init(&B, CC_BINTS, H_IRR, 7, 7, 5, 5, 1, "B <ab|cd>");
    dpd_contract444(&Cmnef, &B, &Sijab, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&B);
    dpd_buf4_close(&Sijab);
    dpd_buf4_init(&X, EOM_TMP, C_irr, 2, 10, 2, 10, 0, "Xijma");
    dpd_buf4_init(&F, CC_FINTS, H_IRR, 10, 7, 10, 5, 1, "F <ia|bc>");
    dpd_contract444(&Cmnef, &F, &X, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&F);
    dpd_buf4_close(&Cmnef);
    dpd_buf4_init(&WM, EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WabefDD_M");
    dpd_file2_init(&tia, CC_OEI, H_IRR, 0, 1, "tia");
    dpd_contract244(&tia, &X, &WM, 0, 2, 1, 1.0, 0.0);
    dpd_file2_close(&tia);
    dpd_buf4_sort(&WM, EOM_TMP, pqsr, 2, 5, "WabefDD_P");
    dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 2, 5, 2, 7, 0, Sijab_lbl);
    dpd_buf4_axpy(&WM, &Sijab, -1.0);
    dpd_buf4_close(&WM);
    dpd_buf4_init(&WP, EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WabefDD_P");
    dpd_buf4_axpy(&WP, &Sijab, 1.0);
    dpd_buf4_close(&WP);
    dpd_buf4_close(&Sijab);
    dpd_buf4_init(&X, EOM_TMP, C_irr, 2, 2, 2, 2, 0, "XIJMN");
    dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 2, 7, 2, 7, 0, Cmnef_lbl);
    dpd_buf4_init(&D, CC_DINTS, H_IRR, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
    dpd_contract444(&Cmnef, &D, &X, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&Cmnef);
    dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 2, 7, 2, 7, 0, Sijab_lbl);
    dpd_buf4_init(&tau, CC_TAMPS, H_IRR, 2, 7, 2, 7, 0, "tauijab");
    dpd_contract444(&X, &tau, &Sijab, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&tau);
    dpd_buf4_close(&Sijab);
    dpd_buf4_close(&X);

#ifdef EOM_DEBUG
  check_sum("WabefDD",i,C_irr);
#endif

    /* SIjAb += WAbEf*CIjEf */
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    sprintf(CMnEf_lbl, "%s %d", "CMnEf", i);
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, CMnEf_lbl);

    /* make use of a more efficient algorithm */
    dpd_buf4_init(&Z, EOM_TMP, C_irr, 5, 0, 5, 0, 0, "Z(Ab,Ij)");
    dpd_buf4_init(&B, CC_BINTS, H_IRR, 5, 5, 5, 5, 0, "B <ab|cd>");
    /*  dpd_contract444(&CMnEf, &B, &SIjAb, 0, 0, 1.0, 1.0); */
    dpd_contract444(&B, &CMnEf, &Z, 0, 0, 1, 0);
    dpd_buf4_close(&B);
    dpd_buf4_sort(&Z, EOM_TMP, rspq, 0, 5, "Z(Ij,Ab)");
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
    dpd_buf4_axpy(&Z, &SIjAb, 1);
    dpd_buf4_close(&Z);

    dpd_buf4_close(&CMnEf);
    dpd_buf4_init(&X, EOM_TMP, C_irr, 0, 10, 0, 10, 0, "XIjMa");
    dpd_buf4_init(&CMnEf, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "CMnfE");
    dpd_buf4_init(&F, CC_FINTS, H_IRR, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_contract444(&CMnEf, &F, &X, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&F);
    dpd_buf4_close(&CMnEf);
    dpd_buf4_sort(&X, EOM_TMP, pqsr, 0, 11, "XIjaM");
    dpd_buf4_close(&X);
    dpd_buf4_init(&X, EOM_TMP, C_irr, 0, 11, 0, 11, 0, "XIjaM");
    dpd_file2_init(&tia, CC_OEI, H_IRR, 0, 1, "tia");
    dpd_contract424(&X, &tia, &SIjAb, 3, 0, 0, -1.0, 1.0);
    dpd_file2_close(&tia);
    dpd_buf4_close(&X);
    dpd_buf4_init(&X, EOM_TMP, C_irr, 0, 10, 0, 10, 0, "XIjMb");
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, CMnEf_lbl);
    dpd_buf4_init(&F, CC_FINTS, H_IRR, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_contract444(&CMnEf, &F, &X, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&F);
    dpd_buf4_close(&CMnEf);
    dpd_file2_init(&tIA, CC_OEI, H_IRR, 0, 1, "tIA");
    dpd_contract244(&tIA, &X, &SIjAb, 0, 2, 1, -1.0, 1.0);
    dpd_file2_close(&tIA);
    dpd_buf4_close(&X);
    dpd_buf4_init(&X, EOM_TMP, C_irr, 0, 0, 0, 0, 0, "XIjMn");
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, CMnEf_lbl);
    dpd_buf4_init(&D, CC_DINTS, H_IRR, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract444(&CMnEf, &D, &X, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&CMnEf);
    dpd_buf4_init(&tau, CC_TAMPS, H_IRR, 0, 5, 0, 5, 0, "tauIjAb");
    dpd_contract444(&X, &tau, &SIjAb, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&tau);
    dpd_buf4_close(&SIjAb);
    dpd_buf4_close(&X);

#ifdef EOM_DEBUG
  check_sum("WabefDD",i,C_irr);
#endif

    sprintf(SIA_lbl, "%s %d", "SIA", i);
    dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, SIA_lbl);
    sprintf(Sia_lbl, "%s %d", "Sia", i);
    dpd_file2_init(&Sia, EOM_Sia, C_irr, 0, 1, Sia_lbl);
    dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 2, 7, 2, 7, 0, SIJAB_lbl);
    dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 2, 7, 2, 7, 0, Sijab_lbl);
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
      c_clean(&SIA,&Sia,&SIJAB,&Sijab,&SIjAb);
    dpd_buf4_close(&SIJAB);
    dpd_buf4_close(&Sijab);
    dpd_buf4_close(&SIjAb);

  }

#ifdef EOM_DEBUG
  check_sum("WabefDD",i,C_irr);
#endif
  return;
}
