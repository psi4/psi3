#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/* This function computes the non-R0 parts of the 2pdm density matrix
   Gciab = 0.5 *(rho_abci + rho_ciab) */

void x_Gciab_6(void);
void x_Gciab_7(void);
void x_Gciab_8(void);
void x_Gciab_9(void);

void x_Gciab(void) { 
  int h, nirreps, i, j, k, a, I, J, K, A, Isym, Jsym, Ksym, Asym, row, col;
  int II,JJ,IIsym,JJsym;
  int L_irr, R_irr, G_irr;
  double value;
  dpdfile2 L1, T1, R1, I1;
  dpdbuf4 G, V, T, L, Z, Z2, R, Tau;
 
  L_irr = params.L_irr; R_irr = params.R_irr; G_irr = params.G_irr;
  nirreps = moinfo.nirreps;

  /* term 1, rho_abci += Lmiab * Rmc */
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 7, 11, 7, 11, 0, "L2R1_VVOV(pqsr)");
  dpd_buf4_sort(&Z, EOM_TMP0, rspq, 11, 7, "GCIAB");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 7, 11, 7, 11, 0, "L2R1_vvov(pqsr)");
  dpd_buf4_sort(&Z, EOM_TMP0, rspq, 11, 7, "Gciab");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 5, 11, 5, 11, 0, "L2R1_VvoV(pqsr)");
  dpd_buf4_sort(&Z, EOM_TMP0, rspq, 11, 5, "GCiAb");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 5, 11, 5, 11, 0, "L2R1_VvOv(pqsr)");
  dpd_buf4_sort(&Z, EOM_TMP0, rsqp, 11, 5, "GcIaB");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 11, 7, 11, 7, 0, "GCIAB");
  dpd_buf4_scm(&Z, -1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 11, 7, 11, 7, 0, "Gciab");
  dpd_buf4_scm(&Z, -1.0);
  dpd_buf4_close(&Z);

  /* term 2, rho_ciab += Rmiab * Lmc */
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 7, 11, 7, 11, 0, "R2L1_VVOV(pqsr)");
  dpd_buf4_sort(&Z, EOM_TMP1, rspq, 11, 7, "GCIAB");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 7, 11, 7, 11, 0, "R2L1_vvov(pqsr)");
  dpd_buf4_sort(&Z, EOM_TMP1, rspq, 11, 7, "Gciab");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 5, 11, 5, 11, 0, "R2L1_VvoV(pqsr)");
  dpd_buf4_sort(&Z, EOM_TMP1, rspq, 11, 5, "GCiAb");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 5, 11, 5, 11, 0, "R2L1_VvOv(pqsr)");
  dpd_buf4_sort(&Z, EOM_TMP1, rsqp, 11, 5, "GcIaB");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 7, 11, 7, 0, "GCIAB");
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 7, 11, 7, 0, "GCIAB");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 7, 11, 7, 0, "Gciab");
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 7, 11, 7, 0, "Gciab");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GCiAb");
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "GCiAb");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GcIaB");
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "GcIaB");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  psio_close(EOM_TMP1,0);
  psio_open(EOM_TMP1, PSIO_OPEN_NEW);

  /* term 3, rho_CIAB -= 0.5 LMNCE RMNAB tIE */
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 11, 2, 11, 2, 0, "Z(CI,MN)");
  dpd_buf4_init(&L, CC_GL, L_irr, 2, 5, 2, 7, 0, "LIJAB");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&L, &T1, &Z, 3, 1, 1, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&L);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 7, 11, 7, 0, "GCIAB");
  dpd_buf4_init(&R, CC_GR, R_irr, 2, 7, 2, 7, 0, "RIJAB");
  dpd_contract444(&Z, &R, &G, 0, 1, -1.0, 1.0);
  dpd_buf4_close(&R);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);
  /* term 3, rho_ciab -= 0.5 Lmnce Rmnab tie */
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 11, 2, 11, 2, 0, "Z(ci,mn)");
  dpd_buf4_init(&L, CC_GL, L_irr, 2, 5, 2, 7, 0, "Lijab");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&L, &T1, &Z, 3, 1, 1, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&L);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 7, 11, 7, 0, "Gciab");
  dpd_buf4_init(&R, CC_GR, R_irr, 2, 7, 2, 7, 0, "Rijab");
  dpd_contract444(&Z, &R, &G, 0, 1, -1.0, 1.0);
  dpd_buf4_close(&R);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);
  /* term 3, rho_CiAb -= 0.5 LMnCe RMnAb tie */
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 11, 0, 11, 0, 0, "Z(Ci,Mn)");
  dpd_buf4_init(&L, CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&L, &T1, &Z, 3, 1, 1, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&L);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GCiAb");
  dpd_buf4_init(&R, CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjAb");
  dpd_contract444(&Z, &R, &G, 0, 1, -1.0, 1.0);
  dpd_buf4_close(&R);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);
  /* term 3, rho_cIaB -= 0.5 LmNcE RmNaB tIE */
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 11, 0, 11, 0, 0, "Z(cI,mN)");
  dpd_buf4_init(&L, CC_GL, L_irr, 0, 5, 0, 5, 0, "LiJaB");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&L, &T1, &Z, 3, 1, 1, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&L);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GcIaB");
  dpd_buf4_init(&R, CC_GR, R_irr, 0, 5, 0, 5, 0, "RiJaB");
  dpd_contract444(&Z, &R, &G, 0, 1, -1.0, 1.0);
  dpd_buf4_close(&R);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);

  psio_close(EOM_TMP1,0);
  psio_open(EOM_TMP1, PSIO_OPEN_NEW);

  /* term 4, rho_CIAB -= 0.5 LMNCE TauMNAB RIE */
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 11, 2, 11, 2, 0, "Z(CI,MN)");
  dpd_buf4_init(&L, CC_GL, L_irr, 2, 5, 2, 7, 0, "LIJAB");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract424(&L, &R1, &Z, 3, 1, 1, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&L);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 7, 11, 7, 0, "GCIAB");
  dpd_buf4_init(&T, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
  dpd_contract444(&Z, &T, &G, 0, 1, -1.0, 1.0);
  dpd_buf4_close(&T);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);
  /* term 4, rho_ciab -= 0.5 Lmnce Taumnab Rie */
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 11, 2, 11, 2, 0, "Z(ci,mn)");
  dpd_buf4_init(&L, CC_GL, L_irr, 2, 5, 2, 7, 0, "Lijab");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
  dpd_contract424(&L, &R1, &Z, 3, 1, 1, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&L);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 7, 11, 7, 0, "Gciab");
  dpd_buf4_init(&T, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
  dpd_contract444(&Z, &T, &G, 0, 1, -1.0, 1.0);
  dpd_buf4_close(&T);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);
  /* term 4, rho_CiAb -= 0.5 LMnCe TauMnAb Rie */
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 11, 0, 11, 0, 0, "Z(Ci,Mn)");
  dpd_buf4_init(&L, CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
  dpd_contract424(&L, &R1, &Z, 3, 1, 1, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&L);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GCiAb");
  dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
  dpd_contract444(&Z, &T, &G, 0, 1, -1.0, 1.0);
  dpd_buf4_close(&T);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);
  /* term 4, rho_cIaB -= 0.5 LmNcE TaumNaB RIE */
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 11, 0, 11, 0, 0, "Z(cI,mN)");
  dpd_buf4_init(&L, CC_GL, L_irr, 0, 5, 0, 5, 0, "LiJaB");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract424(&L, &R1, &Z, 3, 1, 1, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&L);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GcIaB");
  dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauiJaB");
  dpd_contract444(&Z, &T, &G, 0, 1, -1.0, 1.0);
  dpd_buf4_close(&T);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);

  /* term 5, rho_ciab += - (Lmnec Rme) (Tniab + P(ij) Tna Tib) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 7, 11, 7, 0, "GCIAB");
  dpd_buf4_init(&Tau, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tauIJAB");
  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
  dpd_contract244(&I1, &Tau, &G, 0, 0, 0, 1.0, 1.0);
  dpd_file2_close(&I1);
  dpd_buf4_close(&Tau);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 7, 11, 7, 0, "Gciab");
  dpd_buf4_init(&Tau, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tauijab");
  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_ov");
  dpd_contract244(&I1, &Tau, &G, 0, 0, 0, 1.0, 1.0);
  dpd_file2_close(&I1);
  dpd_buf4_close(&Tau);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GCiAb");
  dpd_buf4_init(&Tau, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
  dpd_contract244(&I1, &Tau, &G, 0, 0, 0, 1.0, 1.0);
  dpd_file2_close(&I1);
  dpd_buf4_close(&Tau);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GcIaB");
  dpd_buf4_init(&Tau, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauiJaB");
  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_ov");
  dpd_contract244(&I1, &Tau, &G, 0, 0, 0, 1.0, 1.0);
  dpd_file2_close(&I1);
  dpd_buf4_close(&Tau);
  dpd_buf4_close(&G);

  /* +P(ab) LR_VV(c,a) t(i,b) */
  x_Gciab_6();
  /* +P(ab) LR_TT(c,a) R(i,b) */
  x_Gciab_7();
  /* -P(ab) Lmnce Rinae Tmb */
  x_Gciab_8();
  psio_close(EOM_TMP1,0);
  psio_open(EOM_TMP1, PSIO_OPEN_NEW);
  /* - P(ab) Lmnce Tinae Rmb */
  x_Gciab_9();
  psio_close(EOM_TMP1,0);
  psio_open(EOM_TMP1, PSIO_OPEN_NEW);

  /* term 10, -P(AB) LNMCE TIE TNA RMB */
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 2, 11, 2, 11, 0, "Z(NM,CI)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_buf4_init(&L, CC_GL, L_irr, 2, 5, 2, 7, 0, "LIJAB");
  dpd_contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_buf4_close(&L);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 0, 11, 2, 11, 0, "Z(NM,CI)");
  dpd_buf4_init(&Z2, EOM_TMP1, L_irr, 11, 11, 11, 11, 0, "Z(CI,AM)");
  dpd_contract244(&T1, &Z, &Z2, 0, 0, 1, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "Z(CI,AB)");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract424(&Z2, &R1, &Z, 3, 0, 0, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&Z2);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 7, 0, "GCIAB");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_sort(&Z, EOM_TMP1, pqsr, 11, 5, "Z(CI,BA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "Z(CI,BA)");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  /* term 10, -P(ab) lmnce tie tna rmb */
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 2, 11, 2, 11, 0, "Z(nm,ci)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_buf4_init(&L, CC_GL, L_irr, 2, 5, 2, 7, 0, "Lijab");
  dpd_contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_buf4_close(&L);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 0, 11, 2, 11, 0, "Z(nm,ci)");
  dpd_buf4_init(&Z2, EOM_TMP1, L_irr, 11, 11, 11, 11, 0, "Z(ci,am)");
  dpd_contract244(&T1, &Z, &Z2, 0, 0, 1, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "Z(ci,ab)");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
  dpd_contract424(&Z2, &R1, &Z, 3, 0, 0, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&Z2);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 7, 0, "Gciab");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_sort(&Z, EOM_TMP1, pqsr, 11, 5, "Z(ci,ba)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "Z(ci,ba)");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  /* term 10, GCiAb -= P(AB) LNmCe Tie TNA Rmb */
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 0, 11, 0, 11, 0, "Z(Nm,Ci)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_buf4_init(&L, CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
  dpd_contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_buf4_close(&L);
  dpd_file2_close(&T1);
  dpd_buf4_init(&Z2, EOM_TMP1, L_irr, 11, 11, 11, 11, 0, "Z(Ci,Am)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract244(&T1, &Z, &Z2, 0, 0, 1, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GCiAb");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
  dpd_contract424(&Z2, &R1, &G, 3, 0, 0, -1.0, 1.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&Z2);
  dpd_buf4_close(&G);

  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 0, 11, 0, 11, 0, "Z(nM,Ci)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_buf4_init(&L, CC_GL, L_irr, 0, 5, 0, 5, 0, "LiJAb");
  dpd_contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_buf4_close(&L);
  dpd_file2_close(&T1);
  dpd_buf4_init(&Z2, EOM_TMP1, L_irr, 11, 11, 11, 11, 0, "Z(Ci,bM)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract244(&T1, &Z, &Z2, 0, 0, 1, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "Z(Ci,bA)");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract424(&Z2, &R1, &Z, 3, 0, 0, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&Z2);
  dpd_buf4_sort(&Z, EOM_TMP1, pqsr, 11, 5, "Z(Ci,Ab)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "Z(Ci,Ab)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GCiAb");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  /* term 10, GcIaB - LnMcE TIE Tna RMB + LNmcE TIE TNB Rma */
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 0, 11, 0, 11, 0, "Z(nM,cI)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_buf4_init(&L, CC_GL, L_irr, 0, 5, 0, 5, 0, "LiJaB");
  dpd_contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_buf4_close(&L);
  dpd_file2_close(&T1);
  dpd_buf4_init(&Z2, EOM_TMP1, L_irr, 11, 11, 11, 11, 0, "Z(cI,aM)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract244(&T1, &Z, &Z2, 0, 0, 1, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GcIaB");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract424(&Z2, &R1, &G, 3, 0, 0, -1.0, 1.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&Z2);
  dpd_buf4_close(&G);

  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 0, 11, 0, 11, 0, "Z(Nm,cI)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_buf4_init(&L, CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjaB");
  dpd_contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_buf4_close(&L);
  dpd_file2_close(&T1);
  dpd_buf4_init(&Z2, EOM_TMP1, L_irr, 11, 11, 11, 11, 0, "Z(cI,Bm)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract244(&T1, &Z, &Z2, 0, 0, 1, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "Z(cI,Ba)");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
  dpd_contract424(&Z2, &R1, &Z, 3, 0, 0, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&Z2);
  dpd_buf4_sort(&Z, EOM_TMP1, pqsr, 11, 5, "Z(cI,aB)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "Z(cI,aB)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GcIaB");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);


  /* add 1/2 to ground-state parts of density */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 7, 11, 7, 0, "GCIAB");
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 11, 7, 11, 7, 0, "GCIAB");
  dpd_buf4_axpy(&G, &V, 0.5);
  dpd_buf4_close(&V);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 7, 11, 7, 0, "Gciab");
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 11, 7, 11, 7, 0, "Gciab");
  dpd_buf4_axpy(&G, &V, 0.5);
  dpd_buf4_close(&V);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GCiAb");
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 11, 5, 11, 5, 0, "GCiAb");
  dpd_buf4_axpy(&G, &V, 0.5);
  dpd_buf4_close(&V);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GcIaB");
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 11, 5, 11, 5, 0, "GcIaB");
  dpd_buf4_axpy(&G, &V, 0.5);
  dpd_buf4_close(&V);
  dpd_buf4_close(&G);

  /* clear out temporary files */
  psio_close(EOM_TMP0, 0);
  psio_open(EOM_TMP0, PSIO_OPEN_NEW);
}
