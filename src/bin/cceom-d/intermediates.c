#define EXTERN
#include "globals.h"

/* LR1_OO  =  LIE * RJE */
/* LR1_oo  = Lie * Rje */
/* LR1_VV = LMA * RMB */
/* LR1_vv = Lma * Rmb */
/* L2R1_OV(I,A) = RME * LIMAE + Rme + LImAe */
/* L2R1_OV(i,a) = Rme * Limae + RME + LiMaE */
/* L1R2_OV(I,A) = LME * RIMAE + Lme * RImAe */
/* L1R2_ov(i,a) = Lme * Rimae + LME * RiMaE */
/* LR2_OO(I,J)  = 0.5 * LIMEF * RJMEF + LImEf * RJmEf */
/* LR2_oo(i,j)  = 0.5 * Limef * Rjmef + LiMeF * RjMeF */
/* LR2_VV(A,B) = 0.5 * LMNEA * RMNEB + LmNeA * RmNeB */
/* LR2_vv(a,b) = 0.5 * Lmnea * Rmneb + LMnEa * RMnEb */
/* LT2_OO(I,J) = 0.5 * LIMEF * TJMEF + LImEf * TJmEf */
/* LT2_oo(i,j) = 0.5 * Limef * Tjmef + LiMeF * TjMeF */
/* LT2_VV(A,B) = 0.5 * LMNEA * TMNEB + LmNeA * TmNeB */
/* LT2_vv(a,b) = 0.5 * Lmnea * Tmneb + LMnEa * TMnEb */

/* LR_OO = LR1_OO + LR2_OO */
/* LR_oo = LR1_oo + LR2_oo */
/* LR_VV = LR1_VV + LR2_VV */
/* LR_vv = LR1_vv + LR2_vv */

/* in CC_TAMPS
 tauIJAB  = tIJAB + 1.0 * (TIA * TJB - TIB * TJA)
 tautIJAB = tIJAB + 0.5 * (TIA * TJB - TIB * TJA)
*/

void intermediates(int R_irr)
{
  dpdfile2 L1, R1, T1, I, LR1, LR2;
  dpdbuf4 L2, T2, R2;
  double tval;

  /* LR1_OO(I,J)  =  LIE * RJE */
  dpd_file2_init(&I, EOM_TMP, 0, 0, 0, "LR1_OO");

  dpd_file2_init(&L1, CC_LAMPS, R_irr, 0, 1, "LIA");
  dpd_file2_init(&R1, CC_RAMPS, R_irr, 0, 1, "RIA");
  dpd_contract222(&L1, &R1, &I, 0, 0, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_file2_close(&L1);

  dpd_file2_close(&I);

  /* LR1_oo(i,j)  = Lia * Rje */
  dpd_file2_init(&I, EOM_TMP, 0, 0, 0, "LR1_oo");

  dpd_file2_init(&L1, CC_LAMPS, R_irr, 0, 1, "Lia");
  dpd_file2_init(&R1, CC_RAMPS, R_irr, 0, 1, "Ria");
  dpd_contract222(&L1, &R1, &I, 0, 0, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_file2_close(&L1);

  dpd_file2_close(&I);

  /* LR1_VV(A,B) = LMA * RMB */
  dpd_file2_init(&I, EOM_TMP, 0, 1, 1, "LR1_VV");

  dpd_file2_init(&L1, CC_LAMPS, R_irr, 0, 1, "LIA");
  dpd_file2_init(&R1, CC_RAMPS, R_irr, 0, 1, "RIA");
  dpd_contract222(&L1, &R1, &I, 1, 1, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_file2_close(&L1);

  dpd_file2_close(&I);

  /* LR1_vv(a,b) = Lma * Rmb */
  dpd_file2_init(&I, EOM_TMP, 0, 1, 1, "LR1_vv");

  dpd_file2_init(&L1, CC_LAMPS, R_irr, 0, 1, "Lia");
  dpd_file2_init(&R1, CC_RAMPS, R_irr, 0, 1, "Ria");
  dpd_contract222(&L1, &R1, &I, 1, 1, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_file2_close(&L1);

  dpd_file2_close(&I);

  /* L2R1_OV(I,A) = RME * LIMAE + Rme + LImAe */
  dpd_file2_init(&I, EOM_TMP, 0, 0, 1, "L2R1_OV");

  dpd_buf4_init(&L2, CC_LAMPS, R_irr, 0, 5, 2, 7, 0, "LIJAB");
  dpd_file2_init(&R1, CC_RAMPS, R_irr, 0, 1, "RIA");
  dpd_dot24(&R1, &L2, &I, 0, 0, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&L2);

  dpd_buf4_init(&L2, CC_LAMPS, R_irr, 0, 5, 0, 5, 0, "LIjAb");
  dpd_file2_init(&R1, CC_RAMPS, R_irr, 0, 1, "Ria");
  dpd_dot24(&R1, &L2, &I, 0, 0, 1.0, 1.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&L2);

  dpd_file2_close(&I);

  /* L2R1_OV(i,a) = Rme * Limae + RME + LiMaE */
  dpd_file2_init(&I, EOM_TMP, 0, 0, 1, "L2R1_ov");

  dpd_buf4_init(&L2, CC_LAMPS, R_irr, 0, 5, 2, 7, 0, "Lijab");
  dpd_file2_init(&R1, CC_RAMPS, R_irr, 0, 1, "Ria");
  dpd_dot24(&R1, &L2, &I, 0, 0, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&L2);

  dpd_buf4_init(&L2, CC_LAMPS, R_irr, 0, 5, 0, 5, 0, "LiJaB");
  dpd_file2_init(&R1, CC_RAMPS, R_irr, 0, 1, "RIA");
  dpd_dot24(&R1, &L2, &I, 0, 0, 1.0, 1.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&L2);

  dpd_file2_close(&I);

  /* L1R2_OV(I,A) = LME * RIMAE + Lme * RImAe */
  dpd_file2_init(&I, EOM_TMP, 0, 0, 1, "L1R2_OV");

  dpd_buf4_init(&R2, CC_RAMPS, R_irr, 0, 5, 2, 7, 0, "RIJAB");
  dpd_file2_init(&L1, CC_LAMPS, R_irr, 0, 1, "LIA");
  dpd_dot24(&L1, &R2, &I, 0, 0, 1.0, 0.0);
  dpd_file2_close(&L1);
  dpd_buf4_close(&R2);

  dpd_buf4_init(&R2, CC_RAMPS, R_irr, 0, 5, 0, 5, 0, "RIjAb");
  dpd_file2_init(&L1, CC_LAMPS, R_irr, 0, 1, "Lia");
  dpd_dot24(&L1, &R2, &I, 0, 0, 1.0, 1.0);
  dpd_file2_close(&L1);
  dpd_buf4_close(&R2);

  dpd_file2_close(&I);

  /* L1R2_ov(i,a) = Lme * Rimae + LME * RiMaE */
  dpd_file2_init(&I, EOM_TMP, 0, 0, 1, "L1R2_ov");

  dpd_buf4_init(&R2, CC_RAMPS, R_irr, 0, 5, 2, 7, 0, "Rijab");
  dpd_file2_init(&L1, CC_LAMPS, R_irr, 0, 1, "Lia");
  dpd_dot24(&L1, &R2, &I, 0, 0, 1.0, 0.0);
  dpd_file2_close(&L1);
  dpd_buf4_close(&R2);

  dpd_buf4_init(&R2, CC_RAMPS, R_irr, 0, 5, 0, 5, 0, "RiJaB");
  dpd_file2_init(&L1, CC_LAMPS, R_irr, 0, 1, "LIA");
  dpd_dot24(&L1, &R2, &I, 0, 0, 1.0, 1.0);
  dpd_file2_close(&L1);
  dpd_buf4_close(&R2);

  dpd_file2_close(&I);

  /* LR2_OO(I,J)  = 0.5 * LIMEF * RJMEF + LImEf * RJmEf */
  dpd_file2_init(&I, EOM_TMP, 0, 0, 0, "LR2_OO");

  dpd_buf4_init(&L2, CC_LAMPS, R_irr, 0, 7, 2, 7, 0, "LIJAB");
  dpd_buf4_init(&R2, CC_RAMPS, R_irr, 0, 7, 2, 7, 0, "RIJAB");
  dpd_contract442(&L2, &R2, &I, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&R2);
  dpd_buf4_close(&L2);

  dpd_buf4_init(&L2, CC_LAMPS, R_irr, 0, 5, 0, 5, 0, "LIjAb");
  dpd_buf4_init(&R2, CC_RAMPS, R_irr, 0, 5, 0, 5, 0, "RIjAb");
  dpd_contract442(&L2, &R2, &I, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&R2);
  dpd_buf4_close(&L2);

  dpd_file2_close(&I);

  /* LR2_oo(i,j)  = 0.5 * Limef * Rjmef + LiMeF * RjMeF */
  dpd_file2_init(&I, EOM_TMP, 0, 0, 0, "LR2_oo");

  dpd_buf4_init(&L2, CC_LAMPS, R_irr, 0, 7, 2, 7, 0, "Lijab");
  dpd_buf4_init(&R2, CC_RAMPS, R_irr, 0, 7, 2, 7, 0, "Rijab");
  dpd_contract442(&L2, &R2, &I, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&R2);
  dpd_buf4_close(&L2);

  dpd_buf4_init(&L2, CC_LAMPS, R_irr, 0, 5, 0, 5, 0, "LiJaB");
  dpd_buf4_init(&R2, CC_RAMPS, R_irr, 0, 5, 0, 5, 0, "RiJaB");
  dpd_contract442(&L2, &R2, &I, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&R2);
  dpd_buf4_close(&L2);

  dpd_file2_close(&I);

  /* LR2_VV(A,B) = 0.5 * LMNEA * RMNEB + LmNeA * RmNeB */
  dpd_file2_init(&I, EOM_TMP, 0, 1, 1, "LR2_VV");

  dpd_buf4_init(&L2, CC_LAMPS, R_irr, 2, 5, 2, 7, 0, "LIJAB");
  dpd_buf4_init(&R2, CC_RAMPS, R_irr, 2, 5, 2, 7, 0, "RIJAB");
  dpd_contract442(&L2, &R2, &I, 3, 3, 1.0, 0.0);
  dpd_buf4_close(&R2);
  dpd_buf4_close(&L2);

  dpd_buf4_init(&L2, CC_LAMPS, R_irr, 0, 5, 0, 5, 0, "LiJaB");
  dpd_buf4_init(&R2, CC_RAMPS, R_irr, 0, 5, 0, 5, 0, "RiJaB");
  dpd_contract442(&L2, &R2, &I, 3, 3, 1.0, 1.0);
  dpd_buf4_close(&R2);
  dpd_buf4_close(&L2);

  dpd_file2_close(&I);

  /* LR2_vv(a,b) = 0.5 * Lmnea * Rmneb + LMnEa * RMnEb */
  dpd_file2_init(&I, EOM_TMP, 0, 1, 1, "LR2_vv");

  dpd_buf4_init(&L2, CC_LAMPS, R_irr, 2, 5, 2, 7, 0, "Lijab");
  dpd_buf4_init(&R2, CC_RAMPS, R_irr, 2, 5, 2, 7, 0, "Rijab");
  dpd_contract442(&L2, &R2, &I, 3, 3, 1.0, 0.0);
  dpd_buf4_close(&R2);
  dpd_buf4_close(&L2);

  dpd_buf4_init(&L2, CC_LAMPS, R_irr, 0, 5, 0, 5, 0, "LIjAb");
  dpd_buf4_init(&R2, CC_RAMPS, R_irr, 0, 5, 0, 5, 0, "RIjAb");
  dpd_contract442(&L2, &R2, &I, 3, 3, 1.0, 1.0);
  dpd_buf4_close(&R2);
  dpd_buf4_close(&L2);

  dpd_file2_close(&I);

  /* LT2_OO(I,J) = 0.5 * LIMEF * TJMEF + LImEf * TJmEf */
  dpd_file2_init(&I, EOM_TMP, 0, 0, 0, "LT2_OO");

  dpd_buf4_init(&L2, CC_LAMPS, R_irr, 0, 7, 2, 7, 0, "LIJAB");
  dpd_buf4_init(&T2, CC_TAMPS, R_irr, 0, 7, 2, 7, 0, "tIJAB");
  dpd_contract442(&L2, &T2, &I, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&L2);

  dpd_buf4_init(&L2, CC_LAMPS, R_irr, 0, 5, 0, 5, 0, "LIjAb");
  dpd_buf4_init(&T2, CC_TAMPS, R_irr, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract442(&L2, &T2, &I, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&L2);

  dpd_file2_close(&I);

  /* LT2_oo(i,j) = 0.5 * Limef * Tjmef + LiMeF * TjMeF */
  dpd_file2_init(&I, EOM_TMP, 0, 0, 0, "LT2_oo");

  dpd_buf4_init(&L2, CC_LAMPS, R_irr, 0, 7, 2, 7, 0, "Lijab");
  dpd_buf4_init(&T2, CC_TAMPS, R_irr, 0, 7, 2, 7, 0, "tijab");
  dpd_contract442(&L2, &T2, &I, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&L2);

  dpd_buf4_init(&L2, CC_LAMPS, R_irr, 0, 5, 0, 5, 0, "LiJaB");
  dpd_buf4_init(&T2, CC_TAMPS, R_irr, 0, 5, 0, 5, 0, "tiJaB");
  dpd_contract442(&L2, &T2, &I, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&L2);

  dpd_file2_close(&I);

  /* LT2_VV(A,B) = 0.5 * LMNEA * TMNEB + LmNeA * TmNeB */
  dpd_file2_init(&I, EOM_TMP, 0, 1, 1, "LT2_VV");

  dpd_buf4_init(&L2, CC_LAMPS, R_irr, 2, 5, 2, 7, 0, "LIJAB");
  dpd_buf4_init(&T2, CC_TAMPS, R_irr, 2, 5, 2, 7, 0, "tIJAB");
  dpd_contract442(&L2, &T2, &I, 3, 3, 1.0, 0.0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&L2);

  dpd_buf4_init(&L2, CC_LAMPS, R_irr, 0, 5, 0, 5, 0, "LiJaB");
  dpd_buf4_init(&T2, CC_TAMPS, R_irr, 0, 5, 0, 5, 0, "tiJaB");
  dpd_contract442(&L2, &T2, &I, 3, 3, 1.0, 1.0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&L2);

  dpd_file2_close(&I);

  /* LT2_vv(a,b) = 0.5 * Lmnea * Tmneb + LMnEa * TMnEb */
  dpd_file2_init(&I, EOM_TMP, 0, 1, 1, "LT2_vv");

  dpd_buf4_init(&L2, CC_LAMPS, R_irr, 2, 5, 2, 7, 0, "Lijab");
  dpd_buf4_init(&T2, CC_TAMPS, R_irr, 2, 5, 2, 7, 0, "tijab");
  dpd_contract442(&L2, &T2, &I, 3, 3, 1.0, 0.0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&L2);

  dpd_buf4_init(&L2, CC_LAMPS, R_irr, 0, 5, 0, 5, 0, "LIjAb");
  dpd_buf4_init(&T2, CC_TAMPS, R_irr, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract442(&L2, &T2, &I, 3, 3, 1.0, 1.0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&L2);

  dpd_file2_close(&I);

/* LR_OO = LR1_OO + LR2_OO */
/* LR_oo = LR1_oo + LR2_oo */
/* LR_VV = LR1_VV + LR2_VV */
/* LR_vv = LR1_vv + LR2_vv */
  dpd_file2_init(&I, EOM_TMP, 0, 0, 0, "LR_OO");
  dpd_file2_init(&LR1, EOM_TMP, 0, 0, 0, "LR1_OO");
  dpd_file2_init(&LR2, EOM_TMP, 0, 0, 0, "LR2_OO");
  dpd_file2_axpbycz(&LR1, &LR2, &I, 1.0, 1.0, 0.0);
  dpd_file2_close(&LR2);
  dpd_file2_close(&LR1);
  dpd_file2_close(&I);

  dpd_file2_init(&I, EOM_TMP, 0, 0, 0, "LR_oo");
  dpd_file2_init(&LR1, EOM_TMP, 0, 0, 0, "LR1_oo");
  dpd_file2_init(&LR2, EOM_TMP, 0, 0, 0, "LR2_oo");
  dpd_file2_axpbycz(&LR1, &LR2, &I, 1.0, 1.0, 0.0);
  dpd_file2_close(&LR2);
  dpd_file2_close(&LR1);
  dpd_file2_close(&I);

  dpd_file2_init(&I, EOM_TMP, 0, 1, 1, "LR_VV");
  dpd_file2_init(&LR1, EOM_TMP, 0, 1, 1, "LR1_VV");
  dpd_file2_init(&LR2, EOM_TMP, 0, 1, 1, "LR2_VV");
  dpd_file2_axpbycz(&LR1, &LR2, &I, 1.0, 1.0, 0.0);
  dpd_file2_close(&LR2);
  dpd_file2_close(&LR1);
  dpd_file2_close(&I);

  dpd_file2_init(&I, EOM_TMP, 0, 1, 1, "LR_vv");
  dpd_file2_init(&LR1, EOM_TMP, 0, 1, 1, "LR1_vv");
  dpd_file2_init(&LR2, EOM_TMP, 0, 1, 1, "LR2_vv");
  dpd_file2_axpbycz(&LR1, &LR2, &I, 1.0, 1.0, 0.0);
  dpd_file2_close(&LR2);
  dpd_file2_close(&LR1);
  dpd_file2_close(&I);

  return;
}