#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/* This function computes the non-R0 parts of the 2pdm density matrix
   Gijka = 0.5 *(rho_kaij + rho_Gijka) */

extern void x_Gijka_6(void);
extern void x_Gijka_7(void);
extern void x_Gijka_8(void);
extern void x_Gijka_9(void);

void x_Gijka(void) { 
  int h, nirreps, i, j, k, a, I, J, K, A, Isym, Jsym, Ksym, Asym, row, col;
  int II,JJ,IIsym,JJsym;
  int L_irr, R_irr, G_irr;
  double value;
  dpdfile2 L1A, T1A, L1B, T1B, R1A, R1B, I1A, I1B;
  dpdbuf4 G, V, T, L, Z, Z1, Z2, Tau;
 
  L_irr = params.L_irr; R_irr = params.R_irr; G_irr = params.G_irr;
  nirreps = moinfo.nirreps;

  /* term 1, rho_kaij += Lijae * Rke */
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 2, 10, 2, 10, 0, "L2R1_OOVO(pqsr)");
  dpd_buf4_copy(&Z, EOM_TMP0, "GIJKA");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 2, 10, 2, 10, 0, "L2R1_oovo(pqsr)");
  dpd_buf4_copy(&Z, EOM_TMP0, "Gijka");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 0, 10, 0, 10, 0, "L2R1_OovO(pqsr)");
  dpd_buf4_copy(&Z, EOM_TMP0, "GIjKa");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 0, 10, 0, 10, 0, "L2R1_OoVo(qpsr)");
  dpd_buf4_copy(&Z, EOM_TMP0, "GiJkA");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GIjKa");
  dpd_buf4_scm(&Z, -1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GiJkA");
  dpd_buf4_scm(&Z, -1.0);
  dpd_buf4_close(&Z);

  /* term 2, rho_ijka += Rijae * Lke */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 2, 10, 2, 10, 0, "GIJKA");
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 2, 10, 2, 10, 0, "L1R2_OOVO(pqsr)");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 2, 10, 2, 10, 0, "Gijka");
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 2, 10, 2, 10, 0, "L1R2_oovo(pqsr)");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GIjKa");
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 0, 10, 0, 10, 0, "L1R2_OovO(pqsr)");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GiJkA");
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 0, 10, 0, 10, 0, "L1R2_OoVo(qpsr)");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  
  /* term 3, rho_ijka += 0.5 Rijef Lkmef tma  */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 2, 10, 2, 10, 0, "GIJKA");
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 2, 0, 2, 2, 0, "R2L2_OOOO");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&Z, &T1A, &G, 3, 0, 0, 1.0, 1.0);
  dpd_file2_close(&T1A);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 2, 10, 2, 10, 0, "Gijka");
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 2, 0, 2, 2, 0, "R2L2_oooo");
  dpd_file2_init(&T1B, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&Z, &T1B, &G, 3, 0, 0, 1.0, 1.0);
  dpd_file2_close(&T1B);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GIjKa");
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 0, 0, 0, 0, 0, "R2L2_OoOo");
  dpd_file2_init(&T1B, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&Z, &T1B, &G, 3, 0, 0, 1.0, 1.0);
  dpd_file2_close(&T1B);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GiJkA");
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 0, 0, 0, 0, 0, "R2L2_oOoO");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&Z, &T1A, &G, 3, 0, 0, 1.0, 1.0);
  dpd_file2_close(&T1A);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  /* term 4, rho_ijka += 0.5 [Tijef + P(ij) Tie Tjf] Lkmef Rma */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 2, 10, 2, 10, 0, "GIJKA");
  dpd_buf4_init(&Z, EOM_TMP, L_irr, 2, 0, 2, 2, 0, "Tau2L2_OOOO");
  dpd_file2_init(&R1A, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract424(&Z, &R1A, &G, 3, 0, 0, 1.0, 1.0);
  dpd_file2_close(&R1A);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 2, 10, 2, 10, 0, "Gijka");
  dpd_buf4_init(&Z, EOM_TMP, L_irr, 2, 0, 2, 2, 0, "Tau2L2_oooo");
  dpd_file2_init(&R1B, CC_GR, R_irr, 0, 1, "Ria");
  dpd_contract424(&Z, &R1B, &G, 3, 0, 0, 1.0, 1.0);
  dpd_file2_close(&R1B);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GIjKa");
  dpd_buf4_init(&Z, EOM_TMP, L_irr, 0, 0, 0, 0, 0, "Tau2L2_OoOo");
  dpd_file2_init(&R1B, CC_GR, R_irr, 0, 1, "Ria");
  dpd_contract424(&Z, &R1B, &G, 3, 0, 0, 1.0, 1.0);
  dpd_file2_close(&R1B);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GiJkA");
  dpd_buf4_init(&Z, EOM_TMP, L_irr, 0, 0, 0, 0, 0, "Tau2L2_oOoO");
  dpd_file2_init(&R1A, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract424(&Z, &R1A, &G, 3, 0, 0, 1.0, 1.0);
  dpd_file2_close(&R1A);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  /* term 5, rho_ijka += - (Lkmef Rmf) (Tijea - P(ij) Tie Tja) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 2, 10, 2, 10, 0, "GIJKA");
  dpd_buf4_init(&Tau, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tauIJAB");
  dpd_file2_init(&I1A, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
  dpd_contract244(&I1A, &Tau, &G, 1, 2, 1, -1.0, 1.0);
  dpd_file2_close(&I1A);
  dpd_buf4_close(&Tau);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 2, 10, 2, 10, 0, "Gijka");
  dpd_buf4_init(&Tau, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tauijab");
  dpd_file2_init(&I1B, EOM_TMP, G_irr, 0, 1, "L2R1_ov");
  dpd_contract244(&I1B, &Tau, &G, 1, 2, 1, -1.0, 1.0);
  dpd_file2_close(&I1B);
  dpd_buf4_close(&Tau);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GIjKa");
  dpd_buf4_init(&Tau, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
  dpd_file2_init(&I1A, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
  dpd_contract244(&I1A, &Tau, &G, 1, 2, 1, -1.0, 1.0);
  dpd_file2_close(&I1A);
  dpd_buf4_close(&Tau);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GiJkA");
  dpd_buf4_init(&Tau, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauiJaB");
  dpd_file2_init(&I1B, EOM_TMP, G_irr, 0, 1, "L2R1_ov");
  dpd_contract244(&I1B, &Tau, &G, 1, 2, 1, -1.0, 1.0);
  dpd_file2_close(&I1B);
  dpd_buf4_close(&Tau);
  dpd_buf4_close(&G);

  x_Gijka_6();
  x_Gijka_7();
  /* term 8, +P(ij) Lkmfe rimae tjf */
  x_Gijka_8();
  /* term 9, +P(ij) Lkmfe Timae Rjf, uses Z3, Z4 */
  x_Gijka_9();
  /* term 10, +P(IJ) LKMEF RJF TMA TIE */
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 0, 2, 0, 2, 0, "Z5(JI,KM)");
  dpd_buf4_init(&V, EOM_TMP, G_irr, 2, 10, 2, 10, 0, "L2R1_OOVO(pqsr)");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&V, &T1A, &Z, 3, 1, 1, 1.0, 0.0);
  dpd_buf4_close(&V);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 0, 0, 0, 2, 0, "Z5(JI,KM)");
  dpd_buf4_init(&Z2, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "Z5(JI,KA)");
  dpd_contract424(&Z, &T1A, &Z2, 3, 0, 0, 1.0, 0.0);
  dpd_file2_close(&T1A);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 10, 2, 10, 0, "GIJKA");
  dpd_buf4_axpy(&Z2, &G, -1.0);
  dpd_buf4_sort(&Z2, EOM_TMP0, qprs, 0, 10, "Z5(IJ,KA)");
  dpd_buf4_close(&Z2);
  dpd_buf4_init(&Z2, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "Z5(IJ,KA)");
  dpd_buf4_axpy(&Z2, &G, 1.0);
  dpd_buf4_close(&Z2);
  dpd_buf4_close(&G);
  /* term 10, +P(ij) lkmef rjf tma tie */
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 0, 2, 0, 2, 0, "Z5(ji,km)");
  dpd_buf4_init(&V, EOM_TMP, G_irr, 2, 10, 2, 10, 0, "L2R1_oovo(pqsr)");
  dpd_file2_init(&T1B, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&V, &T1B, &Z, 3, 1, 1, 1.0, 0.0);
  dpd_buf4_close(&V);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 0, 0, 0, 2, 0, "Z5(ji,km)");
  dpd_buf4_init(&Z2, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "Z5(ji,ka)");
  dpd_contract424(&Z, &T1B, &Z2, 3, 0, 0, 1.0, 0.0);
  dpd_file2_close(&T1B);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 10, 2, 10, 0, "Gijka");
  dpd_buf4_axpy(&Z2, &G, -1.0);
  dpd_buf4_sort(&Z2, EOM_TMP0, qprs, 0, 10, "Z5(ij,ka)");
  dpd_buf4_close(&Z2);
  dpd_buf4_init(&Z2, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "Z5(ij,ka)");
  dpd_buf4_axpy(&Z2, &G, 1.0);
  dpd_buf4_close(&Z2);
  dpd_buf4_close(&G);

  /* term 10, GIjKa += LKmEf Rjf TIE tma */
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 0, 0, 0, 0, 0, "Z5(Ij,Km)");
  dpd_buf4_init(&V, EOM_TMP, G_irr, 0, 11, 0, 11, 0, "L2R1_OoVo");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract244(&T1A, &V, &Z, 1, 2, 0, 1.0, 0.0);
  dpd_file2_close(&T1A);
  dpd_buf4_close(&V);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GIjKa");
  dpd_file2_init(&T1B, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&Z, &T1B, &G, 3, 0, 0, 1.0, 1.0);
  dpd_file2_close(&T1B);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 0, 0, 0, 0, 0, "Z5(Ij,Km)");
  dpd_buf4_init(&V, EOM_TMP, G_irr, 0, 10, 0, 10, 0, "L2R1_OovO(pqsr)");
  dpd_file2_init(&T1B, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&V, &T1B, &Z, 3, 1, 1, 1.0, 0.0);
  dpd_file2_close(&T1B);
  dpd_buf4_close(&V);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GIjKa");
  dpd_file2_init(&T1B, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&Z, &T1B, &G, 3, 0, 0, 1.0, 1.0);
  dpd_file2_close(&T1B);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  /* term 10, GiJkA += P(ij) LkMeF RJF Tie tMA */
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 0, 0, 0, 0, 0, "Z5(iJ,Km)");
  dpd_buf4_init(&V, EOM_TMP, G_irr, 0, 11, 0, 11, 0, "L2R1_OovO(qprs)");
  dpd_file2_init(&T1B, CC_OEI, 0, 0, 1, "tia");
  dpd_contract244(&T1B, &V, &Z, 1, 2, 0, 1.0, 0.0);
  dpd_file2_close(&T1B);
  dpd_buf4_close(&V);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GiJkA");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&Z, &T1A, &G, 3, 0, 0, 1.0, 1.0);
  dpd_file2_close(&T1A);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 0, 0, 0, 0, 0, "Z6(iJ,kM)");
  dpd_buf4_init(&V, EOM_TMP, G_irr, 0, 10, 0, 10, 0, "L2R1_OoVo(qpsr)");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&V, &T1A, &Z, 3, 1, 1, 1.0, 0.0);
  dpd_file2_close(&T1A);
  dpd_buf4_close(&V);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GiJkA");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&Z, &T1A, &G, 3, 0, 0, 1.0, 1.0);
  dpd_file2_close(&T1A);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  /* add 1/2 to ground-state parts of density */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 2, 10, 2, 10, 0, "GIJKA");
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 2, 10, 2, 10, 0, "GIJKA");
  dpd_buf4_axpy(&G, &V, 0.5);
  dpd_buf4_close(&V);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 2, 10, 2, 10, 0, "Gijka");
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 2, 10, 2, 10, 0, "Gijka");
  dpd_buf4_axpy(&G, &V, 0.5);
  dpd_buf4_close(&V);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GIjKa");
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 0, 10, 0, 10, 0, "GIjKa");
  dpd_buf4_axpy(&G, &V, 0.5);
  dpd_buf4_close(&V);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GiJkA");
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 0, 10, 0, 10, 0, "GiJkA");
  dpd_buf4_axpy(&G, &V, 0.5);
  dpd_buf4_close(&V);
  dpd_buf4_close(&G);

  /* clear out temporary files */
  psio_close(EOM_TMP0, 0);
  psio_open(EOM_TMP0, PSIO_OPEN_NEW);
}
