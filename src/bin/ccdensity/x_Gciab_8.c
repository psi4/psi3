#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/* This function computes terms of excited Gciab
   term 8  +P(ab) Lmnce Rinae Tmb
   term 9, +P(AB) LMNCE TINAE RMB
*/

void x_Gciab_8(void) { 
  int h, nirreps, i, j, k, a, I, J, K, A, Isym, Jsym, Ksym, Asym, row, col;
  int II,JJ,IIsym,JJsym;
  int L_irr, R_irr, G_irr;
  double value;
  dpdfile2 L1A, T1A, L1B, T1B, R1A, R1B, I1A, I1B;
  dpdbuf4 G, V, T, L, Z, Z1, Z2, Tau;
 
  L_irr = params.L_irr; R_irr = params.R_irr; G_irr = params.G_irr;
  nirreps = moinfo.nirreps;

  /* term 8, +P(AB) LMNCE RINAE TMB */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 5, 10, 5, 0, "Z(IA,BC)");
  dpd_buf4_init(&V, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVOV");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract244(&T1A, &V, &Z, 0, 2, 1, 1.0, 0.0); 
  dpd_file2_close(&T1A);
  dpd_buf4_close(&V);
  dpd_buf4_sort(&Z, EOM_TMP1, sprq, 11, 5, "Z(CI,BA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "Z(CI,BA)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 7, 0, "GCIAB");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_sort(&Z, EOM_TMP1, pqsr, 11, 5, "Z(CI,AB)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "Z(CI,AB)");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);
  /* term 8, +P(ab) Lmnce Rinae Tmb */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 5, 10, 5, 0, "Z(ia,bc)");
  dpd_buf4_init(&V, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_ovov");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tia");
  dpd_contract244(&T1A, &V, &Z, 0, 2, 1, 1.0, 0.0); 
  dpd_file2_close(&T1A);
  dpd_buf4_close(&V);
  dpd_buf4_sort(&Z, EOM_TMP1, sprq, 11, 5, "Z(ci,ba)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "Z(ci,ba)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 7, 0, "Gciab");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_sort(&Z, EOM_TMP1, pqsr, 11, 5, "Z(ci,ab)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "Z(ci,ab)");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);
  
  /* term 8, GCiAb -= LmNCe RiNAe tmb */
  dpd_buf4_init(&V, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_oVoV");
  dpd_buf4_sort(&V, EOM_TMP1, spqr, 11, 11, "Z(Ci,Am)");
  dpd_buf4_close(&V);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 11, 11, 11, 0, "Z(Ci,Am)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GCiAb");
  dpd_file2_init(&T1B, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&Z, &T1B, &G, 3, 0, 0, -1.0, 1.0); 
  dpd_file2_close(&T1B);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);
  /* term 8, GCiAb -= (LMnCe Rinbe + LMNCE RiNbE) TMA */
  dpd_buf4_init(&V, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_ovOV");
  dpd_buf4_sort(&V, EOM_TMP1, sprq, 11, 10, "Z(Ci,Mb)");
  dpd_buf4_close(&V);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 10, 11, 10, 0, "Z(Ci,Mb)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GCiAb");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract244(&T1A, &Z, &G, 0, 2, 1, 1.0, 1.0); 
  dpd_file2_close(&T1A);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);

  /* term 8, GcIaB -= LMncE RInaE tMB */
  dpd_buf4_init(&V, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OvOv");
  dpd_buf4_sort(&V, EOM_TMP1, spqr, 11, 11, "Z(cI,aM)");
  dpd_buf4_close(&V);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 11, 11, 11, 0, "Z(cI,aM)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GcIaB");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&Z, &T1A, &G, 3, 0, 0, -1.0, 1.0); 
  dpd_file2_close(&T1A);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);
  /* term 8, GcIaB -= (LmNcE RINBE + Lmnce RInBe) Tma */
  dpd_buf4_init(&V, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVov");
  dpd_buf4_sort(&V, EOM_TMP1, sprq, 11, 10, "Z(cI,mB)");
  dpd_buf4_close(&V);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 10, 11, 10, 0, "Z(cI,mB)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GcIaB");
  dpd_file2_init(&T1B, CC_OEI, 0, 0, 1, "tia");
  dpd_contract244(&T1B, &Z, &G, 0, 2, 1, 1.0, 1.0); 
  dpd_file2_close(&T1B);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);

  psio_close(EOM_TMP1, 0);
  psio_open(EOM_TMP1, PSIO_OPEN_NEW);

  /* term 9, +P(AB) LMNCE TINAE RMB */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 5, 10, 5, 0, "Z(IA,BC)");
  dpd_buf4_init(&V, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "VIAJB");
  dpd_file2_init(&R1A, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract244(&R1A, &V, &Z, 0, 2, 1, 1.0, 0.0); 
  dpd_file2_close(&R1A);
  dpd_buf4_close(&V);
  dpd_buf4_sort(&Z, EOM_TMP1, sprq, 11, 5, "Z(CI,BA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "Z(CI,BA)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 7, 0, "GCIAB");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_sort(&Z, EOM_TMP1, pqsr, 11, 5, "Z(CI,AB)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "Z(CI,AB)");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);
  /* term 9, +P(ab) Lmnce Tinae Rmb */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 5, 10, 5, 0, "Z(ia,bc)");
  dpd_buf4_init(&V, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "Viajb");
  dpd_file2_init(&R1B, CC_GR, R_irr, 0, 1, "Ria");
  dpd_contract244(&R1B, &V, &Z, 0, 2, 1, 1.0, 0.0); 
  dpd_file2_close(&R1B);
  dpd_buf4_close(&V);
  dpd_buf4_sort(&Z, EOM_TMP1, sprq, 11, 5, "Z(ci,ba)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "Z(ci,ba)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 7, 0, "Gciab");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_sort(&Z, EOM_TMP1, pqsr, 11, 5, "Z(ci,ab)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "Z(ci,ab)");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);
  
  /* term 9, GCiAb -= LmNCe TiNAe Rmb */
  dpd_buf4_init(&V, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "ViAjB");
  dpd_buf4_sort(&V, EOM_TMP1, spqr, 11, 11, "Z(Ci,Am)");
  dpd_buf4_close(&V);
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 11, 11, 11, 11, 0, "Z(Ci,Am)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GCiAb");
  dpd_file2_init(&R1B, CC_GR, R_irr, 0, 1, "Ria");
  dpd_contract424(&Z, &R1B, &G, 3, 0, 0, -1.0, 1.0); 
  dpd_file2_close(&R1B);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);
  /* term 9, GCiAb -= (LMnCe Tinbe + LMNCE TiNbE) RMA */
  dpd_buf4_init(&V, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "ViaJB");
  dpd_buf4_sort(&V, EOM_TMP1, sprq, 11, 10, "Z(Ci,Mb)");
  dpd_buf4_close(&V);
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 11, 10, 11, 10, 0, "Z(Ci,Mb)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GCiAb");
  dpd_file2_init(&R1A, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract244(&R1A, &Z, &G, 0, 2, 1, 1.0, 1.0); 
  dpd_file2_close(&R1A);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);

  /* term 9, GcIaB -= LMncE TInaE RMB */
  dpd_buf4_init(&V, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "VIaJb");
  dpd_buf4_sort(&V, EOM_TMP1, spqr, 11, 11, "Z(cI,aM)");
  dpd_buf4_close(&V);
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 11, 11, 11, 11, 0, "Z(cI,aM)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GcIaB");
  dpd_file2_init(&R1A, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract424(&Z, &R1A, &G, 3, 0, 0, -1.0, 1.0); 
  dpd_file2_close(&R1A);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);
  /* term 9, GcIaB -= (LmNcE TINBE + Lmnce TInBe) Rma */
  dpd_buf4_init(&V, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "VIAjb");
  dpd_buf4_sort(&V, EOM_TMP1, sprq, 11, 10, "Z(cI,mB)");
  dpd_buf4_close(&V);
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 11, 10, 11, 10, 0, "Z(cI,mB)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GcIaB");
  dpd_file2_init(&R1B, CC_GR, R_irr, 0, 1, "Ria");
  dpd_contract244(&R1B, &Z, &G, 0, 2, 1, 1.0, 1.0); 
  dpd_file2_close(&R1B);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);

  return;
}
