#include "globals.h"

/* This function computes the extra contributions to sigma_1 and sigma_2
  for EOM_CC3 computations that are not normally present in a EOM_CCSD
  calculation */

/* The additional terms are
 * <S| H    <T| (Uhat C2)c   |0> |T> / (w-wt) -> sigma_1
 * <D| H    <T| (Uhat C2)c   |0> |T> / (w-wt) -> sigma_2
 *
 * <S| Hhat <T| (Utilde T2)c |0> |T> / (w-wt) -> sigma_1
 * <D| Hhat <T| (Utilde T2)c |0> |T> / (w-wt) -> sigma_2
 *
 * <D| H'   <T| (Uhat T2)c   |0> |T> / (-wt) -> sigma_2
 *
 *  See Eqn. (83) in JCP, 103, 7429, 1995
 *
 *  All three terms can be evaluated by the same functions given
 *  different matrix elements.
 *
 * */

extern void cc3_sigma_UHF_AAA(int term, dpdbuf4 *CMNEF, dpdfile2 *SIA, dpdbuf4 *SIJAB, double omega);
extern void cc3_sigma_UHF_BBB(int term, dpdbuf4 *Cmnef, dpdfile2 *Sia, dpdbuf4 *Sijab, double omega);

void cc3_sigma_UHF_AAB(int term, dpdbuf4 *CMNEF, dpdbuf4 *CMnEf, dpdbuf4 *CmNeF,
        dpdfile2 *SIA, dpdfile2 *Sia, dpdbuf4 *SIJAB, dpdbuf4 *SIjAb);

extern void cc3_sigma_UHF_BBA(int term);

void sigmaCC3(int i, int C_irr, double omega) {
  dpdfile2 CME, Cme, SIA, Sia;
  dpdbuf4 CMNEF, CMnEf, Cmnef, SIJAB, Sijab, SIjAb;
  dpdbuf4 TIJAB;
  char lbl[32];




  if (params.eom_ref == 0) { /* RHF */
  }
  else if (params.eom_ref == 1) { /* ROHF */
  }
  else if (params.eom_ref == 2) { /* UHF */

    /*
      * <S| H    <T| (Uhat C2)c   |0> |T> / (w-wt) -> sigma_1
      * <D| Hhat <T| (Uhat C2)c   |0> |T> / (w-wt) -> sigma_2
    */

    sprintf(lbl, "%s %d", "SIA", i);
    dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);

    sprintf(lbl, "%s %d", "SIJAB", i);
    dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 0, 5, 2, 7, 0, lbl);

    sprintf(lbl, "%s %d", "CMNEF", i);
    dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 0, 5, 2, 7, 0, lbl);

    cc3_sigma_UHF_AAA(1, &CMNEF, &SIA, &SIJAB, 0.0);

    dpd_buf4_close(&CMNEF);

    /*
      * <S| H    <T| (Uhat C2)c   |0> |T> / (w-wt) -> sigma_1
      * <D| Hhat <T| (Uhat C2)c   |0> |T> / (w-wt) -> sigma_2
    */
    dpd_buf4_init(&TIJAB, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");

    cc3_sigma_UHF_AAA(2, &TIJAB, &SIA, &SIJAB, 0.0);

    dpd_buf4_close(&TIJAB);


    /*
     * * <D| H'   <T| (Uhat T2)c   |0> |T> / (-wt) -> sigma_2
    */

    cc3_sigma_UHF_AAA(3, &TIJAB, &SIA, &SIJAB, 0.0);


    dpd_file2_close(&SIA);
    dpd_buf4_close(&SIJAB);

  }

  return;
}


