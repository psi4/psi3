#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/* This function computes 
   term 8 of Gijka */

void x_Gijka_8(void) { 
  int h, nirreps, i, j, k, a, I, J, K, A, Isym, Jsym, Ksym, Asym, row, col;
  int II,JJ,IIsym,JJsym;
  int L_irr, R_irr, G_irr;
  double value;
  dpdfile2 L1A, T1A, L1B, T1B, R1A, R1B, I1A, I1B;
  dpdbuf4 G, V, T, L, Z, Z1, Z2, Tau;
 
  L_irr = params.L_irr; R_irr = params.R_irr; G_irr = params.G_irr;
  nirreps = moinfo.nirreps;

  /* term 8, +P(ij) Lkmfe rimae tjf */
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 10, 0, 10, 0, 0, "Z(IA,KJ)");
  dpd_buf4_init(&V, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVOV");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&V, &T1A, &Z, 3, 1, 0, 1.0, 0.0); 
  dpd_file2_close(&T1A);
  dpd_buf4_close(&V);
  dpd_buf4_sort(&Z, EOM_TMP0, psrq, 0, 10, "Z(IJ,KA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "Z(IJ,KA)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 10, 2, 10, 0, "GIJKA");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_sort(&Z, EOM_TMP0, qprs, 0, 10, "Z(JI,KA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "Z(JI,KA)");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);

  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 10, 0, 10, 0, 0, "Z(ia,kj)");
  dpd_buf4_init(&V, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_ovov");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&V, &T1A, &Z, 3, 1, 0, 1.0, 0.0); 
  dpd_file2_close(&T1A);
  dpd_buf4_close(&V);
  dpd_buf4_sort(&Z, EOM_TMP0, psrq, 0, 10, "Z(ij,ka)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "Z(ij,ka)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 10, 2, 10, 0, "Gijka");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_sort(&Z, EOM_TMP0, qprs, 0, 10, "Z(ji,ka)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "Z(ji,ka)");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);

  /* GIjKa += R2L2_OvOv(Ia,Kf) T(j,f) */
  /* GIjKa -= R2L2_OvOv(ja,KF) T(I,F) */
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 10, 0, 10, 0, 0, "Z(Ia,Kj)");
  dpd_buf4_init(&V, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OvOv");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&V, &T1A, &Z, 3, 1, 0, 1.0, 0.0); 
  dpd_file2_close(&T1A);
  dpd_buf4_close(&V);
  dpd_buf4_sort(&Z, EOM_TMP0, psrq, 0, 10, "Z(Ij,Ka)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "Z(Ij,Ka)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GIjKa");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 10, 0, 10, 0, 0, "Z2(ja,KI)");
  dpd_buf4_init(&V, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_ovOV");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&V, &T1A, &Z, 3, 1, 0, 1.0, 0.0); 
  dpd_file2_close(&T1A);
  dpd_buf4_close(&V);
  dpd_buf4_sort(&Z, EOM_TMP0, psrq, 0, 10, "Z2(jI,Ka)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "Z2(jI,Ka)");
  dpd_buf4_sort(&Z, EOM_TMP0, qprs, 0, 10, "Z2(Ij,Ka)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "Z2(Ij,Ka)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GIjKa");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  /* GiJkA += R2L2_oVoV(iA,kF) T(J,F) */
  /* GiJkA += R2L2_OVov(JA,kf) T(i,f) */
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 10, 0, 10, 0, 0, "Z(iA,kJ)");
  dpd_buf4_init(&V, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_oVoV");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&V, &T1A, &Z, 3, 1, 0, 1.0, 0.0); 
  dpd_file2_close(&T1A);
  dpd_buf4_close(&V);
  dpd_buf4_sort(&Z, EOM_TMP0, psrq, 0, 10, "Z(iJ,kA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "Z(iJ,kA)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GiJkA");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 10, 0, 10, 0, 0, "Z2(JA,ki)");
  dpd_buf4_init(&V, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVov");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&V, &T1A, &Z, 3, 1, 0, 1.0, 0.0); 
  dpd_file2_close(&T1A);
  dpd_buf4_close(&V);
  dpd_buf4_sort(&Z, EOM_TMP0, psrq, 0, 10, "Z2(Ji,kA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "Z2(Ji,kA)");
  dpd_buf4_sort(&Z, EOM_TMP0, qprs, 0, 10, "Z2(iJ,kA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "Z2(iJ,kA)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GiJkA");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  return;
}
