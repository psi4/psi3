#define EXTERN
#include "globals.h"

/* build one-electron intermediates for the construction of the xi amplitudes */

void x_xi_intermediates(void)
{
  dpdfile2 L1, R1, T1, I, LR1, LR2, LT1, LT2;
  dpdbuf4 L2, T2, R2, D, H2, I2;
  int L_irr, R_irr, G_irr;

  L_irr = params.L_irr;
  R_irr = params.R_irr;
  G_irr = params.G_irr;

  /* RD_OO(I,J)  =  RIMEF * <JM||EF> + RImEf <Jm|Ef> */
  dpd_file2_init(&I, EOM_TMP_XI, R_irr, 0, 0, "RD_OO");
  dpd_buf4_init(&R2, CC_GR, R_irr, 0, 7, 2, 7, 0, "RIJAB");
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
  dpd_contract442(&R2, &D, &I, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&D);
  dpd_buf4_close(&R2);
  dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjAb");
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_contract442(&R2, &D, &I, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&D);
  dpd_buf4_close(&R2);
  dpd_file2_close(&I);

  /* RD_oo(i,j)  =  Rimef * <jm||ef> + RiMeF <jM|eF> */
  dpd_file2_init(&I, EOM_TMP_XI, R_irr, 0, 0, "RD_oo");
  dpd_buf4_init(&R2, CC_GR, R_irr, 0, 7, 2, 7, 0, "Rijab");
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
  dpd_contract442(&R2, &D, &I, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&D);
  dpd_buf4_close(&R2);
  dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 0, 5, 0, "RiJaB");
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_contract442(&R2, &D, &I, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&D);
  dpd_buf4_close(&R2);
  dpd_file2_close(&I);

  /* RD_VV(E,A)  =  RMNFE * <MN||FA> + RmNfE <mN|fA> */
  dpd_file2_init(&I, EOM_TMP_XI, R_irr, 1, 1, "RD_VV");
  dpd_buf4_init(&R2, CC_GR, R_irr, 2, 5, 2, 7, 0, "RIJAB");
  dpd_buf4_init(&D, CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");
  dpd_contract442(&R2, &D, &I, 3, 3, 1.0, 0.0);
  dpd_buf4_close(&D);
  dpd_buf4_close(&R2);
  dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 0, 5, 0, "RiJaB");
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_contract442(&R2, &D, &I, 3, 3, 1.0, 1.0);
  dpd_buf4_close(&D);
  dpd_buf4_close(&R2);
  dpd_file2_close(&I);

  /* RD_vv(e,a)  =  Rmnfe * <mn||fa> + RMnFe <Mn|Fa> */
  dpd_file2_init(&I, EOM_TMP_XI, R_irr, 1, 1, "RD_vv");
  dpd_buf4_init(&R2, CC_GR, R_irr, 2, 5, 2, 7, 0, "Rijab");
  dpd_buf4_init(&D, CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");
  dpd_contract442(&R2, &D, &I, 3, 3, 1.0, 0.0);
  dpd_buf4_close(&D);
  dpd_buf4_close(&R2);
  dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjAb");
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_contract442(&R2, &D, &I, 3, 3, 1.0, 1.0);
  dpd_buf4_close(&D);
  dpd_buf4_close(&R2);
  dpd_file2_close(&I);

  /* R2Wamef_OOVO = Rmnfg Weifg -> OOOV (mn,ie) */
  dpd_buf4_init(&I2, EOM_TMP, R_irr, 2, 10, 2, 10, 0, "R2Wamef_OOOV");
  dpd_buf4_init(&R2, CC_GR, R_irr, 2, 7, 2, 7, 0, "RIJAB");
  dpd_buf4_init(&H2, CC_HBAR, 0, 10, 7, 10, 7, 0, "WAMEF");
  dpd_contract444(&R2, &H2, &I2, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&H2);
  dpd_buf4_close(&R2);
  dpd_buf4_close(&I2);

  dpd_buf4_init(&I2, EOM_TMP, R_irr, 0, 10, 0, 10, 0, "R2Wamef_OoOv");
  dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjaB");
  dpd_buf4_init(&H2, CC_HBAR, 0, 10, 5, 10, 5, 0, "WaMeF");
  dpd_contract444(&R2, &H2, &I2, 0, 0, -1.0, 0.0);
  dpd_buf4_close(&H2);
  dpd_buf4_close(&R2);
  dpd_buf4_close(&I2);

  dpd_buf4_init(&I2, EOM_TMP, R_irr, 0, 10, 0, 10, 0, "R2Wamef_oOoV");
  dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 0, 5, 0, "RiJAb");
  dpd_buf4_init(&H2, CC_HBAR, 0, 10, 5, 10, 5, 0, "WAmEf");
  dpd_contract444(&R2, &H2, &I2, 0, 0, -1.0, 0.0);
  dpd_buf4_close(&H2);
  dpd_buf4_close(&R2);
  dpd_buf4_close(&I2);

  dpd_buf4_init(&I2, EOM_TMP, R_irr, 2, 10, 2, 10, 0, "R2Wamef_ooov");
  dpd_buf4_init(&R2, CC_GR, R_irr, 2, 7, 2, 7, 0, "Rijab");
  dpd_buf4_init(&H2, CC_HBAR, 0, 10, 7, 10, 7, 0, "Wamef");
  dpd_contract444(&R2, &H2, &I2, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&H2);
  dpd_buf4_close(&R2);
  dpd_buf4_close(&I2);

  /* R2Wamef_OV(n,g) = Rnmef Wgmef */
  dpd_file2_init(&I, EOM_TMP, R_irr, 0, 1, "R2Wamef_OV");

  dpd_buf4_init(&R2, CC_GR, R_irr, 0, 7, 2, 7, 0, "RIJAB");
  dpd_buf4_init(&H2, CC_HBAR, 0, 11, 7, 11, 7, 0, "WAMEF (AM,E>F)");
  dpd_contract442(&R2, &H2, &I, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&H2);
  dpd_buf4_close(&R2);
  dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjAb");
  dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf (Am,Ef)");
  dpd_contract442(&R2, &H2, &I, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&H2);
  dpd_buf4_close(&R2);

  dpd_file2_close(&I);

  dpd_file2_init(&I, EOM_TMP, R_irr, 0, 1, "R2Wamef_ov");

  dpd_buf4_init(&R2, CC_GR, R_irr, 0, 7, 2, 7, 0, "Rijab");
  dpd_buf4_init(&H2, CC_HBAR, 0, 11, 7, 11, 7, 0, "Wamef (am,e>f)");
  dpd_contract442(&R2, &H2, &I, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&H2);
  dpd_buf4_close(&R2);
  dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 0, 5, 0, "RiJaB");
  dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 5, 0, "WaMeF (aM,eF)");
  dpd_contract442(&R2, &H2, &I, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&H2);
  dpd_buf4_close(&R2);

  dpd_file2_close(&I);

  /* R1Wmnie_OO (M,I) = Rne Wmnie */
  dpd_file2_init(&I, EOM_TMP_XI, R_irr, 0, 0, "R1Wmnie_OO");

  dpd_buf4_init(&H2, CC_HBAR, 0, 0, 10, 2, 10, 0, "WMNIE (M>N,IE)");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
  dpd_dot24(&R1, &H2, &I, 0, 0, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&H2);
  dpd_buf4_init(&H2, CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe (Mn,Ie)");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
  dpd_dot24(&R1, &H2, &I, 0, 0, 1.0, 1.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&H2);

  dpd_file2_close(&I);
  dpd_file2_init(&I, EOM_TMP_XI, R_irr, 0, 0, "R1Wmnie_oo");

  dpd_buf4_init(&H2, CC_HBAR, 0, 0, 10, 2, 10, 0, "Wmnie (m>n,ie)");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
  dpd_dot24(&R1, &H2, &I, 0, 0, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&H2);
  dpd_buf4_init(&H2, CC_HBAR, 0, 0, 10, 0, 10, 0, "WmNiE (mN,iE)");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
  dpd_dot24(&R1, &H2, &I, 0, 0, 1.0, 1.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&H2);

  dpd_file2_close(&I);

  /* R1Wamef_VV (a,e) = Rmf Wamef */
  dpd_file2_init(&I, EOM_TMP_XI, R_irr, 1, 1, "R1Wamef_VV");

  dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 7, 0, "WAMEF (AM,E>F)");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
  dpd_dot24(&R1, &H2, &I, 0, 0, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&H2);
  dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf (Am,Ef)");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
  dpd_dot24(&R1, &H2, &I, 0, 0, 1.0, 1.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&H2);

  dpd_file2_close(&I);
  dpd_file2_init(&I, EOM_TMP_XI, R_irr, 1, 1, "R1Wamef_vv");

  dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 7, 0, "Wamef (am,e>f)");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
  dpd_dot24(&R1, &H2, &I, 0, 0, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&H2);
  dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 5, 0, "WaMeF (aM,eF)");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
  dpd_dot24(&R1, &H2, &I, 0, 0, 1.0, 1.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&H2);

  dpd_file2_close(&I);

  return;
}
