#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/* This function computes term 9 of Gijka,  +P(ij) Lkmfe Timae Rjf */

void x_Gijka_9(void) { 
  int nirreps, L_irr, R_irr, G_irr;
  dpdfile2 R1A, R1B;
  dpdbuf4 G, V, Z;
  double L0;
 
  L_irr = params.L_irr; R_irr = params.R_irr; G_irr = params.G_irr;
  nirreps = moinfo.nirreps; L0 = params.L0;

  /* term 9, +P(ij) Lkmfe rimae tjf */
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 10, 0, 10, 0, 0, "Z3(IA,KJ)");
  dpd_buf4_init(&V, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "VIAJB");
  dpd_file2_init(&R1A, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract424(&V, &R1A, &Z, 3, 1, 0, 1.0, 0.0); 
  dpd_file2_close(&R1A);
  dpd_buf4_close(&V);
  dpd_buf4_sort(&Z, EOM_TMP0, psrq, 0, 10, "Z3(IJ,KA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "Z3(IJ,KA)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 10, 2, 10, 0, "GIJKA");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_sort(&Z, EOM_TMP0, qprs, 0, 10, "Z3(JI,KA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "Z3(JI,KA)");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);

  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 10, 0, 10, 0, 0, "Z3(ia,kj)");
  dpd_buf4_init(&V, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "Viajb");
  dpd_file2_init(&R1B, CC_GR, R_irr, 0, 1, "Ria");
  dpd_contract424(&V, &R1B, &Z, 3, 1, 0, 1.0, 0.0); 
  dpd_file2_close(&R1B);
  dpd_buf4_close(&V);
  dpd_buf4_sort(&Z, EOM_TMP0, psrq, 0, 10, "Z3(ij,ka)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "Z3(ij,ka)");
  dpd_buf4_init(&G, EOM_TMP0, L_irr, 0, 10, 2, 10, 0, "Gijka");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_sort(&Z, EOM_TMP0, qprs, 0, 10, "Z3(ji,ka)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "Z3(ji,ka)");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);

  /* GIjKa += R2L2_OvOv(Ia,Kf) T(j,f) */
  /* GIjKa -= R2L2_OvOv(ja,KF) T(I,F) */
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 10, 0, 10, 0, 0, "Z3(Ia,Kj)");
  dpd_buf4_init(&V, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "VIaJb");
  dpd_file2_init(&R1B, CC_GR, R_irr, 0, 1, "Ria");
  dpd_contract424(&V, &R1B, &Z, 3, 1, 0, 1.0, 0.0); 
  dpd_file2_close(&R1B);
  dpd_buf4_close(&V);
  dpd_buf4_sort(&Z, EOM_TMP0, psrq, 0, 10, "Z3(Ij,Ka)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "Z3(Ij,Ka)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GIjKa");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 10, 0, 10, 0, 0, "Z4(ja,KI)");
  dpd_buf4_init(&V, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "ViaJB");
  dpd_file2_init(&R1A, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract424(&V, &R1A, &Z, 3, 1, 0, 1.0, 0.0); 
  dpd_file2_close(&R1A);
  dpd_buf4_close(&V);
  dpd_buf4_sort(&Z, EOM_TMP0, psrq, 0, 10, "Z4(jI,Ka)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "Z4(jI,Ka)");
  dpd_buf4_sort(&Z, EOM_TMP0, qprs, 0, 10, "Z4(Ij,Ka)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "Z4(Ij,Ka)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GIjKa");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  /* GiJkA += R2L2_oVoV(iA,kF) T(J,F) */
  /* GiJkA += R2L2_OVov(JA,kf) T(i,f) */
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 10, 0, 10, 0, 0, "Z3(iA,kJ)");
  dpd_buf4_init(&V, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "ViAjB");
  dpd_file2_init(&R1A, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract424(&V, &R1A, &Z, 3, 1, 0, 1.0, 0.0); 
  dpd_file2_close(&R1A);
  dpd_buf4_close(&V);
  dpd_buf4_sort(&Z, EOM_TMP0, psrq, 0, 10, "Z3(iJ,kA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "Z3(iJ,kA)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GiJkA");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 10, 0, 10, 0, 0, "Z4(JA,ki)");
  dpd_buf4_init(&V, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "VIAjb");
  dpd_file2_init(&R1B, CC_GR, R_irr, 0, 1, "Ria");
  dpd_contract424(&V, &R1B, &Z, 3, 1, 0, 1.0, 0.0); 
  dpd_file2_close(&R1B);
  dpd_buf4_close(&V);
  dpd_buf4_sort(&Z, EOM_TMP0, psrq, 0, 10, "Z4(Ji,kA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "Z4(Ji,kA)");
  dpd_buf4_sort(&Z, EOM_TMP0, qprs, 0, 10, "Z4(iJ,kA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "Z4(iJ,kA)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GiJkA");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  return;
}
