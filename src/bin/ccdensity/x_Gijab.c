#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

extern void x_Gijab_ROHF_2(void);

void x_Gijab_ROHF(void)
{
  int h, nirreps;
  int R_irr, L_irr, G_irr;
  double value, tval;
  dpdfile2 T1, L1, I1, T1A, T1B, Z1, R1;
  dpdbuf4 R, I, G, L, T, V, Z, Z2;

  nirreps = moinfo.nirreps;
  R_irr = params.R_irr; L_irr = params.L_irr; G_irr = params.G_irr;

  /* (1-R0) * Tau(IJ,AB), if G is totally symmetric like tau */
  if (G_irr == 0) {
    dpd_buf4_init(&T, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    //dpd_buf4_scmcopy(&T, EOM_TMP0, "GIJAB", 1.0 - params.R0);
//    if (params.use_zeta)
//      dpd_buf4_scmcopy(&T, EOM_TMP0, "GIJAB", 1.0+params.RZ_overlap);
//    else
      dpd_buf4_scmcopy(&T, EOM_TMP0, "GIJAB", 1.0);
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
    //dpd_buf4_scmcopy(&T, EOM_TMP0, "Gijab", 1.0 - params.R0);
//    if (params.use_zeta)
//      dpd_buf4_scmcopy(&T, EOM_TMP0, "Gijab", 1.0+params.RZ_overlap);
//    else
      dpd_buf4_scmcopy(&T, EOM_TMP0, "Gijab", 1.0);
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    //dpd_buf4_scmcopy(&T, EOM_TMP0, "GIjAb", 1.0 - params.R0);
//    if (params.use_zeta)
//      dpd_buf4_scmcopy(&T, EOM_TMP0, "GIjAb", 1.0+params.RZ_overlap);
//    else
      dpd_buf4_scmcopy(&T, EOM_TMP0, "GIjAb", 1.0);
    dpd_buf4_close(&T);
  }
  else {
    dpd_buf4_init(&G, EOM_TMP0, G_irr, 2, 7, 2, 7, 0, "GIJAB");
    dpd_buf4_scm(&G,0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, EOM_TMP0, G_irr, 2, 7, 2, 7, 0, "Gijab");
    dpd_buf4_scm(&G,0.0);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 5, 0, 5, 0, "GIjAb");
    dpd_buf4_scm(&G,0.0);
    dpd_buf4_close(&G);
  }

  /* -P(ij) LR_OO(M,I) Tau(MJ,AB); terms 4,5,8,9 */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 7, 0, 7, 0, "Z(IJ,A>B)");
  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 0, "LR_OO");
  // -P(ij) L2R1_OV(M,F) Tau(MJ,AB); terms 35, 36 */
  dpd_file2_init(&Z1, EOM_TMP1, G_irr, 0, 0, "Z(N,I)");
  dpd_file2_axpy(&I1, &Z1, 1.0, 0);
  dpd_file2_close(&I1);
  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract222(&I1, &T1, &Z1, 0, 0, 1.0, 1.0);
  dpd_file2_close(&T1);
  dpd_file2_close(&I1);
  dpd_buf4_init(&T, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tauIJAB");
  dpd_contract244(&Z1, &T, &Z, 0, 0, 0, 1.0, 0.0);
  dpd_file2_close(&Z1);
  dpd_buf4_close(&T);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 7, 2, 7, 0, "GIJAB");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_sort(&Z, EOM_TMP1, qprs, 0, 7, "Z(JI,A>B)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 7, 0, 7, 0, "Z(JI,A>B)");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  /* -P(ij) LR_oo(m,i) Tau(mj,ab); terms 4,5,8,9 */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 7, 0, 7, 0, "Z(ij,a>b)");
  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 0, "LR_oo");
  /* -P(ij) L2R1_ov(m,f) Tau(mj,ab); terms 35, 36 */
  dpd_file2_init(&Z1, EOM_TMP1, G_irr, 0, 0, "Z(n,i)");
  dpd_file2_axpy(&I1, &Z1, 1.0, 0);
  dpd_file2_close(&I1);
  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_ov");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract222(&I1, &T1, &Z1, 0, 0, 1.0, 1.0);
  dpd_file2_close(&T1);
  dpd_file2_close(&I1);
  dpd_buf4_init(&T, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tauijab");
  dpd_contract244(&Z1, &T, &Z, 0, 0, 0, 1.0, 0.0);
  dpd_file2_close(&Z1);
  dpd_buf4_close(&T);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 7, 2, 7, 0, "Gijab");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_sort(&Z, EOM_TMP1, qprs, 0, 7, "Z(ji,a>b)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 7, 0, 7, 0, "Z(ji,a>b)");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  psio_close(EOM_TMP1, 0);
  psio_open(EOM_TMP1,PSIO_OPEN_NEW);
  /* GIjAb += -P(Ij) LR_OO(M,I) Tau(Mj,Ab); terms 4,5,8,9 */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 5, 0, 5, 0, "GIjAb");
  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 0, "LR_OO");
  /* -P(Ij) L2R1_OV(M,F) Tau(Nm,Ab); terms 35, 36 */
  dpd_file2_init(&Z1, EOM_TMP1, G_irr, 0, 0, "Z(N,I)");
  dpd_file2_axpy(&I1, &Z1, 1.0, 0);
  dpd_file2_close(&I1);
  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract222(&I1, &T1, &Z1, 0, 0, 1.0, 1.0);
  dpd_file2_close(&T1);
  dpd_file2_close(&I1);
  dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
  dpd_contract244(&Z1, &T, &G, 0, 0, 0, -1.0, 1.0);
  dpd_file2_close(&Z1);
  dpd_buf4_close(&T);

  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z(jI,bA)");
  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 0, "LR_oo");
  dpd_file2_init(&Z1, EOM_TMP1, G_irr, 0, 0, "Z(n,i)");
  dpd_file2_axpy(&I1, &Z1, 1.0, 0);
  dpd_file2_close(&I1);
  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_ov");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract222(&I1, &T1, &Z1, 0, 0, 1.0, 1.0);
  dpd_file2_close(&T1);
  dpd_file2_close(&I1);
  dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauiJaB");
  dpd_contract244(&Z1, &T, &Z, 0, 0, 0, 1.0, 0.0);
  dpd_file2_close(&Z1);
  dpd_buf4_close(&T);
  dpd_buf4_sort(&Z, EOM_TMP1, qpsr, 0, 5, "Z(Ij,Ab)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&G);

  psio_close(EOM_TMP1,0);
  psio_open(EOM_TMP1,PSIO_OPEN_NEW);

  /* -P(ab) LR_VV(F,A) Tau(IJ,FB); terms 6,7,10,11 */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 2, 5, 2, 5, 0, "Z(I>J,AB)");
  dpd_buf4_init(&T, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tauIJAB");
  dpd_file2_init(&I1, EOM_TMP, G_irr, 1, 1, "LR_VV");
  dpd_contract244(&I1, &T, &Z, 0, 2, 1, 1.0, 0.0);
  dpd_file2_close(&I1);
  dpd_buf4_close(&T);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 2, 5, 2, 7, 0, "GIJAB");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_sort(&Z, EOM_TMP1, pqsr, 2, 5, "Z(I>J,BA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 2, 5, 2, 5, 0, "Z(I>J,BA)");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  /* -P(ab) LR_vv(f,a) Tau(ij,fb); */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 2, 5, 2, 5, 0, "Z(i>j,ab)");
  dpd_buf4_init(&T, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tauijab");
  dpd_file2_init(&I1, EOM_TMP, G_irr, 1, 1, "LR_vv");
  dpd_contract244(&I1, &T, &Z, 0, 2, 1, 1.0, 0.0);
  dpd_file2_close(&I1);
  dpd_buf4_close(&T);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 2, 5, 2, 7, 0, "Gijab");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_sort(&Z, EOM_TMP1, pqsr, 2, 5, "Z(i>j,ba)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 2, 5, 2, 5, 0, "Z(i>j,ba)");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  /* GIjAb += - LR_VV(F,A) Tau(Ij,Fb) + LR_VV(f,a) Tau(Ij,fA) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 5, 0, 5, 0, "GIjAb");
  dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
  dpd_file2_init(&I1, EOM_TMP, G_irr, 1, 1, "LR_VV");
  dpd_contract244(&I1, &T, &G, 0, 2, 1, -1.0, 1.0);
  dpd_file2_close(&I1);
  dpd_buf4_close(&T);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z(Ij,bA)");
  dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjbA");
  dpd_file2_init(&I1, EOM_TMP, G_irr, 1, 1, "LR_vv");
  dpd_contract244(&I1, &T, &Z, 0, 2, 1, 1.0, 0.0);
  dpd_file2_close(&I1);
  dpd_buf4_close(&T);
  dpd_buf4_sort(&Z, EOM_TMP1, pqsr, 0, 5, "Z(Ij,Ab)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&G);

  psio_close(EOM_TMP1,0);
  psio_open(EOM_TMP1,PSIO_OPEN_NEW);

  /* + 1/4 Lmnef Rmnab Tau_ijef, terms 13, 15 */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 2, 7, 2, 7, 0, "GIJAB");
  dpd_buf4_init(&I, EOM_TMP, L_irr, 2, 2, 2, 2, 0, "Tau2L2_OOOO");
  dpd_buf4_init(&R, CC_GR, R_irr, 2, 7, 2, 7, 0, "RIJAB");
  dpd_contract444(&I, &R, &G, 0, 1, 1.0, 1.0);
  dpd_buf4_close(&R);
  dpd_buf4_close(&I);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 2, 7, 2, 7, 0, "Gijab");
  dpd_buf4_init(&I, EOM_TMP, L_irr, 2, 2, 2, 2, 0, "Tau2L2_oooo");
  dpd_buf4_init(&R, CC_GR, R_irr, 2, 7, 2, 7, 0, "Rijab");
  dpd_contract444(&I, &R, &G, 0, 1, 1.0, 1.0);
  dpd_buf4_close(&R);
  dpd_buf4_close(&I);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 5, 0, 5, 0, "GIjAb");
  dpd_buf4_init(&I, EOM_TMP, L_irr, 0, 0, 0, 0, 0, "Tau2L2_OoOo");
  dpd_buf4_init(&R, CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjAb");
  dpd_contract444(&I, &R, &G, 0, 1, 1.0, 1.0);
  dpd_buf4_close(&R);
  dpd_buf4_close(&I);
  dpd_buf4_close(&G);

  /* - 0.5 P(ij) (Lmnfe Rie) Tjf (taumnab), terms 24, 26 */
  /* + 1/4 Lmnef Rijef Tau_mnab, terms 12, 14 */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 2, 0, 2, 0, "Z(IJ,M>N)");
  dpd_buf4_init(&I, EOM_TMP, G_irr, 2, 10, 2, 10, 0, "L2R1_OOVO(pqsr)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&I, &T1, &Z, 3, 1, 1, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&I);
  /* add terms 12, 14 */
  dpd_buf4_init(&I, EOM_TMP, G_irr, 0, 2, 2, 2, 0, "R2L2_OOOO");
  dpd_buf4_axpy(&I, &Z, -0.5);
  dpd_buf4_close(&I);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 0, 7, 0, 7, 0, "Z(IJ,A>B)");
  dpd_buf4_init(&T, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
  dpd_contract444(&Z, &T, &Z2, 0, 1, 1.0, 0.0);
  dpd_buf4_close(&T);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 7, 2, 7, 0, "GIJAB");
  dpd_buf4_axpy(&Z2, &G, -1.0);
  dpd_buf4_sort(&Z2, EOM_TMP1, qprs, 0, 7, "Z(JI,A>B)");
  dpd_buf4_close(&Z2);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 0, 7, 0, 7, 0, "Z(JI,A>B)");
  dpd_buf4_axpy(&Z2, &G, 1.0);
  dpd_buf4_close(&Z2);
  dpd_buf4_close(&G);

  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 2, 0, 2, 0, "Z(ij,m>n)");
  dpd_buf4_init(&I, EOM_TMP, G_irr, 2, 10, 2, 10, 0, "L2R1_oovo(pqsr)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&I, &T1, &Z, 3, 1, 1, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&I);
  /* add terms 12, 14 */
  dpd_buf4_init(&I, EOM_TMP, G_irr, 0, 2, 2, 2, 0, "R2L2_oooo");
  dpd_buf4_axpy(&I, &Z, -0.5);
  dpd_buf4_close(&I);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 0, 7, 0, 7, 0, "Z(ij,a>b)");
  dpd_buf4_init(&T, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
  dpd_contract444(&Z, &T, &Z2, 0, 1, 1.0, 0.0);
  dpd_buf4_close(&T);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 7, 2, 7, 0, "Gijab");
  dpd_buf4_axpy(&Z2, &G, -1.0);
  dpd_buf4_sort(&Z2, EOM_TMP1, qprs, 0, 7, "Z(ji,a>b)");
  dpd_buf4_close(&Z2);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 0, 7, 0, 7, 0, "Z(ji,a>b)");
  dpd_buf4_axpy(&Z2, &G, 1.0);
  dpd_buf4_close(&Z2);
  dpd_buf4_close(&G);

  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 0, 0, 0, 0, "Z(Ij,Mn)");
  dpd_buf4_init(&I, EOM_TMP, G_irr, 0, 10, 0, 10, 0, "L2R1_OovO(pqsr)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&I, &T1, &Z, 3, 1, 1, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&I);
  dpd_buf4_init(&I, EOM_TMP, G_irr, 0, 11, 0, 11, 0, "L2R1_OoVo");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract244(&T1, &I, &Z, 1, 2, 0, 1.0, 1.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&I);
  /* add terms 12, 14 */
  dpd_buf4_init(&I, EOM_TMP, G_irr, 0, 0, 0, 0, 0, "R2L2_OoOo");
  dpd_buf4_axpy(&I, &Z, 1.0);
  dpd_buf4_close(&I);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 5, 0, 5, 0, "GIjAb");
  dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
  dpd_contract444(&Z, &T, &G, 0, 1, 1.0, 1.0);
  dpd_buf4_close(&T);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  psio_close(EOM_TMP1,0);
  psio_open(EOM_TMP1,PSIO_OPEN_NEW);

  /* + 0.5 P(AB) (tau_IJEF LMNEF) RMA TNB ; terms 25, 27 */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 2, 11, 2, 11, 0, "Z(I>J,AN)");
  dpd_buf4_init(&I, EOM_TMP, L_irr, 2, 0, 2, 2, 0, "Tau2L2_OOOO");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract244(&R1, &I, &Z, 0, 2, 1, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&I);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 2, 5, 2, 5, 0, "Z(I>J,AB)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&Z, &T1, &Z2, 3, 0, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 2, 5, 2, 7, 0, "GIJAB");
  dpd_buf4_axpy(&Z2, &G, 1.0);
  dpd_buf4_sort(&Z2, EOM_TMP1, pqsr, 2, 5, "Z(I>J,BA)");
  dpd_buf4_close(&Z2);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 2, 5, 2, 5, 0, "Z(I>J,BA)");
  dpd_buf4_axpy(&Z2, &G, -1.0);
  dpd_buf4_close(&Z2);
  dpd_buf4_close(&G);
  /* + 0.5 P(ab) (tau_ijef Lmnef) Rma Tnb ; terms 25, 27 */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 2, 11, 2, 11, 0, "Z(i>j,an)");
  dpd_buf4_init(&I, EOM_TMP, L_irr, 2, 0, 2, 2, 0, "Tau2L2_oooo");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
  dpd_contract244(&R1, &I, &Z, 0, 2, 1, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&I);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 2, 5, 2, 5, 0, "Z(i>j,ab)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&Z, &T1, &Z2, 3, 0, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 2, 5, 2, 7, 0, "Gijab");
  dpd_buf4_axpy(&Z2, &G, 1.0);
  dpd_buf4_sort(&Z2, EOM_TMP1, pqsr, 2, 5, "Z(i>j,ba)");
  dpd_buf4_close(&Z2);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 2, 5, 2, 5, 0, "Z(i>j,ba)");
  dpd_buf4_axpy(&Z2, &G, -1.0);
  dpd_buf4_close(&Z2);
  dpd_buf4_close(&G);
  /* + tau_IjEf LMnEf RMA Tnb ; terms 25, 27 */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 11, 0, 11, 0, "Z(Ij,An)");
  dpd_buf4_init(&I, EOM_TMP, L_irr, 0, 0, 0, 0, 0, "Tau2L2_OoOo");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract244(&R1, &I, &Z, 0, 2, 1, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&I);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 5, 0, 5, 0, "GIjAb");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&Z, &T1, &G, 3, 0, 0, 1.0, 1.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  /* + tau_IjEf LNmEf Rmb TNA ; terms 25, 27 */
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 0, 11, 0, 11, 0, "Z2(Ij,An)");
  dpd_buf4_init(&I, EOM_TMP, L_irr, 0, 0, 0, 0, 0, "Tau2L2_OoOo");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract244(&T1, &I, &Z, 0, 2, 1, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&I);
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
  dpd_contract424(&Z, &R1, &G, 3, 0, 0, 1.0, 1.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&G);

  psio_close(EOM_TMP1,0);
  psio_open(EOM_TMP1,PSIO_OPEN_NEW);

  /* terms combined to P(ij)P(ab) Z(i,a) T(j,b), 3,22,23,33 */
  x_Gijab_ROHF_2();
  /* terms combined to P(ij)P(ab) Z(i,a) R(j,b), 18,32,34,19 */
  x_Gijab_ROHF_3();

  /* -P(ij)(Lme Tie + 0.5 Lmnef Tinef) Rmjab, term 16, 30 */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 7, 0, 7, 0, "Z(IJ,A>B)");
  dpd_buf4_init(&R, CC_GR, R_irr, 0, 7, 2, 7, 0, "RIJAB");
  dpd_file2_init(&I1, EOM_TMP, L_irr, 0, 0, "LT_OO");
  dpd_contract244(&I1, &R, &Z, 0, 0, 0, 1.0, 0.0);
  dpd_file2_close(&I1);
  dpd_buf4_close(&R);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 7, 2, 7, 0, "GIJAB");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_sort(&Z, EOM_TMP1, qprs, 0, 7, "Z(JI,A>B)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 7, 0, 7, 0, "Z(JI,A>B)");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 7, 0, 7, 0, "Z(ij,a>b)");
  dpd_buf4_init(&R, CC_GR, R_irr, 0, 7, 2, 7, 0, "Rijab");
  dpd_file2_init(&I1, EOM_TMP, L_irr, 0, 0, "LT_oo");
  dpd_contract244(&I1, &R, &Z, 0, 0, 0, 1.0, 0.0);
  dpd_file2_close(&I1);
  dpd_buf4_close(&R);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 7, 2, 7, 0, "Gijab");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_sort(&Z, EOM_TMP1, qprs, 0, 7, "Z(ji,a>b)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 7, 0, 7, 0, "Z(ji,a>b)");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 5, 0, 5, 0, "GIjAb");
  dpd_buf4_init(&R, CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjAb");
  dpd_file2_init(&I1, EOM_TMP, L_irr, 0, 0, "LT_OO");
  dpd_contract244(&I1, &R, &G, 0, 0, 0, -1.0, 1.0);
  dpd_file2_close(&I1);
  dpd_buf4_close(&R);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z(jI,Ab)");
  dpd_buf4_init(&R, CC_GR, R_irr, 0, 5, 0, 5, 0, "RiJAb");
  dpd_file2_init(&I1, EOM_TMP, L_irr, 0, 0, "LT_oo");
  dpd_contract244(&I1, &R, &Z, 0, 0, 0, 1.0, 0.0);
  dpd_file2_close(&I1);
  dpd_buf4_close(&R);
  dpd_buf4_sort(&Z, EOM_TMP1, qprs, 0, 5, "Z(Ij,Ab)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  psio_close(EOM_TMP1,0);
  psio_open(EOM_TMP1, PSIO_OPEN_NEW);

  /* -P(ab)(Lme Tmb + 0.5 Lmnfe Tmnfb) Rijae, term 17, 31 */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 2, 5, 2, 5, 0, "Z(I>J,AB)");
  dpd_buf4_init(&R, CC_GR, R_irr, 2, 5, 2, 7, 0, "RIJAB");
  dpd_file2_init(&I1, EOM_TMP, L_irr, 1, 1, "LT_VV");
  dpd_contract424(&R, &I1, &Z, 3, 0, 0, 1.0, 0.0);
  dpd_file2_close(&I1);
  dpd_buf4_close(&R);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 2, 5, 2, 7, 0, "GIJAB");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_sort(&Z, EOM_TMP1, pqsr, 2, 5, "Z(I>J,BA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 2, 5, 2, 5, 0, "Z(I>J,BA)");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 2, 5, 2, 5, 0, "Z(i>j,ab)");
  dpd_buf4_init(&R, CC_GR, R_irr, 2, 5, 2, 7, 0, "Rijab");
  dpd_file2_init(&I1, EOM_TMP, L_irr, 1, 1, "LT_vv");
  dpd_contract424(&R, &I1, &Z, 3, 0, 0, 1.0, 0.0);
  dpd_file2_close(&I1);
  dpd_buf4_close(&R);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 2, 5, 2, 7, 0, "Gijab");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_sort(&Z, EOM_TMP1, pqsr, 2, 5, "Z(i>j,ba)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 2, 5, 2, 5, 0, "Z(i>j,ba)");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 5, 0, 5, 0, "GIjAb");
  dpd_buf4_init(&R, CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjAb");
  dpd_file2_init(&I1, EOM_TMP, L_irr, 1, 1, "LT_vv");
  dpd_contract424(&R, &I1, &G, 3, 0, 0, -1.0, 1.0);
  dpd_file2_close(&I1);
  dpd_buf4_close(&R);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z(Ij,bA)");
  dpd_buf4_init(&R, CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjaB");
  dpd_file2_init(&I1, EOM_TMP, L_irr, 1, 1, "LT_VV");
  dpd_contract424(&R, &I1, &Z, 3, 0, 0, 1.0, 0.0);
  dpd_file2_close(&I1);
  dpd_buf4_close(&R);
  dpd_buf4_sort(&Z, EOM_TMP1, pqsr, 0, 5, "Z(Ij,Ab)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  psio_close(EOM_TMP1,0);
  psio_open(EOM_TMP1, PSIO_OPEN_NEW);

  /* -P(ab) lmnef rme tijfb tna, term 37 */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 2, 11, 2, 11, 0, "Z(I>J,AN)");
  dpd_buf4_init(&T, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
  dpd_contract424(&T, &I1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&I1);
  dpd_buf4_close(&T);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 2, 5, 2, 5, 0, "Z(I>J,AB)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&Z, &T1, &Z2, 3, 0, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&Z2);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 2, 5, 2, 7, 0, "GIJAB");
  dpd_buf4_init(&Z, EOM_TMP1, 0, 2, 5, 2, 5, 0, "Z(I>J,AB)");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_sort(&Z, EOM_TMP1, pqsr, 2, 5, "Z(I>J,BA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, 0, 2, 5, 2, 5, 0, "Z(I>J,BA)");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_close(&G);

  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 2, 11, 2, 11, 0, "Z(i>j,an)");
  dpd_buf4_init(&T, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");
  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_ov");
  dpd_contract424(&T, &I1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&I1);
  dpd_buf4_close(&T);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 2, 5, 2, 5, 0, "Z(i>j,ab)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&Z, &T1, &Z2, 3, 0, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&Z2);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 2, 5, 2, 7, 0, "Gijab");
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 2, 5, 2, 5, 0, "Z(i>j,ab)");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_sort(&Z, EOM_TMP1, pqsr, 2, 5, "Z(i>j,ba)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 2, 5, 2, 5, 0, "Z(i>j,ba)");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_close(&G);

  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 11, 0, 11, 0, "Z(Ij,An)");
  dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_ov");
  dpd_contract424(&T, &I1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&I1);
  dpd_buf4_close(&T);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 5, 0, 5, 0, "GIjAb");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&Z, &T1, &G, 3, 0, 0, -1.0, 1.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 10, 0, 10, 0, "Z(Ij,Nb)");
  dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
  dpd_contract244(&I1, &T, &Z, 1, 2, 1, 1.0, 0.0);
  dpd_file2_close(&I1);
  dpd_buf4_close(&T);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 5, 0, 5, 0, "GIjAb");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract244(&T1, &Z, &G, 0, 2, 1, -1.0, 1.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  /* compute Z(IA,JB) Z(ia,jb) and Z(IA,jb) for terms
     20, 28, 29, 21 then permute and add in */
  /* + P(ij) P(ab) (Rimae Lnmfe) Tnjfb, term 20 */
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "Z(IA,JB)");
  dpd_buf4_init(&I, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVOV");
  dpd_buf4_init(&T, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
  dpd_contract444(&I, &T, &Z, 0, 1, 1.0, 0.0);
  dpd_buf4_close(&I);
  dpd_buf4_close(&T);
  dpd_buf4_init(&I, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVov");
  dpd_buf4_init(&T, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
  dpd_contract444(&I, &T, &Z, 0, 1, 1.0, 1.0);
  dpd_buf4_close(&I);
  dpd_buf4_close(&T);

  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "Z(ia,jb)");
  dpd_buf4_init(&I, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_ovov");
  dpd_buf4_init(&T, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
  dpd_contract444(&I, &T, &Z, 0, 1, 1.0, 0.0);
  dpd_buf4_close(&I);
  dpd_buf4_close(&T);
  dpd_buf4_init(&I, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_ovOV");
  dpd_buf4_init(&T, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
  dpd_contract444(&I, &T, &Z, 0, 1, 1.0, 1.0);
  dpd_buf4_close(&I);
  dpd_buf4_close(&T);

  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "Z(IA,jb)");
  dpd_buf4_init(&I, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVOV");
  dpd_buf4_init(&T, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
  dpd_contract444(&I, &T, &Z, 0, 1, 1.0, 0.0);
  dpd_buf4_close(&I);
  dpd_buf4_close(&T);
  dpd_buf4_init(&I, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVov");
  dpd_buf4_init(&T, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
  dpd_contract444(&I, &T, &Z, 0, 1, 1.0, 1.0);
  dpd_buf4_close(&I);
  dpd_buf4_close(&T);
  dpd_buf4_init(&I, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_ovOV");
  dpd_buf4_init(&T, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
  dpd_contract444(&T, &I, &Z, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&I);
  dpd_buf4_close(&T);
  dpd_buf4_init(&I, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_ovov");
  dpd_buf4_init(&T, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
  dpd_contract444(&T, &I, &Z, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&I);
  dpd_buf4_close(&T);
  dpd_buf4_close(&Z);

  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z(Ib,jA)");
  dpd_buf4_init(&I, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_oVoV");
  dpd_buf4_init(&T, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
  dpd_contract444(&T, &I, &Z, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&I);
  dpd_buf4_close(&T);
  dpd_buf4_init(&I, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OvOv");
  dpd_buf4_init(&T, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
  dpd_contract444(&I, &T, &Z, 0, 1, 1.0, 1.0);
  dpd_buf4_close(&I);
  dpd_buf4_close(&T);
  dpd_buf4_sort_axpy(&Z, EOM_TMP0, psrq, 10, 10, "Z(IA,jb)", 1.0);
  dpd_buf4_close(&Z);

  psio_close(EOM_TMP1,0);
  psio_open(EOM_TMP1,PSIO_OPEN_NEW);

  /* - P(ij) P(ab) (Tjmbe Lnmfe) Tif Rna, term 28 */
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 10, 0, 10, 0, 0, "Z(JB,NI)");
  dpd_buf4_init(&I, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "VIAJB");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&I, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&I);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 11, 10, 11, 10, 0, "Z(AI,JB)");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract244(&R1, &Z, &Z2, 0, 2, 0, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&Z);
  dpd_buf4_sort_axpy(&Z2, EOM_TMP0, qprs, 10, 10, "Z(IA,JB)", -1.0);
  dpd_buf4_close(&Z2);
  
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 10, 0, 10, 0, 0, "Z(jb,ni)");
  dpd_buf4_init(&I, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "Viajb");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&I, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&I);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 11, 10, 11, 10, 0, "Z(ai,jb)");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
  dpd_contract244(&R1, &Z, &Z2, 0, 2, 0, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&Z);
  dpd_buf4_sort_axpy(&Z2, EOM_TMP0, qprs, 10, 10, "Z(ia,jb)", -1.0);
  dpd_buf4_close(&Z2);

  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 10, 0, 10, 0, 0, "Z(jb,NI)");
  dpd_buf4_init(&I, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "ViaJB");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&I, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&I);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 11, 10, 11, 10, 0, "Z(AI,jb)");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract244(&R1, &Z, &Z2, 0, 2, 0, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&Z);
  dpd_buf4_sort_axpy(&Z2, EOM_TMP0, qprs, 10, 10, "Z(IA,jb)", -1.0);
  dpd_buf4_close(&Z2);

  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 10, 0, 10, 0, 0, "Z(Ib,Nj)");
  dpd_buf4_init(&I, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "VIaJb");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&I, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&I);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(Ib,Aj)");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract244(&R1, &Z, &Z2, 0, 2, 1, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&Z);
  dpd_buf4_sort_axpy(&Z2, EOM_TMP0, prsq, 10, 10, "Z(IA,jb)", +1.0);
  dpd_buf4_close(&Z2);

  psio_close(EOM_TMP1,0);
  psio_open(EOM_TMP1,PSIO_OPEN_NEW);

  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 10, 0, 10, 0, 0, "Z(jA,nI)");
  dpd_buf4_init(&I, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "ViAjB");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&I, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&I);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 11, 10, 11, 10, 0, "Z(bI,jA)");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
  dpd_contract244(&R1, &Z, &Z2, 0, 2, 0, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&Z);
  dpd_buf4_sort(&Z2, EOM_TMP1, sqrp, 11, 10, "Z(AI,jb)");
  dpd_buf4_close(&Z2);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 11, 10, 11, 10, 0, "Z(AI,jb)");
  dpd_buf4_sort_axpy(&Z2, EOM_TMP0, qprs, 10, 10, "Z(IA,jb)", +1.0);
  dpd_buf4_close(&Z2);

  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 10, 0, 10, 0, 0, "Z(IA,nj)");
  dpd_buf4_init(&I, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "VIAjb");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&I, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&I);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(IA,bj)");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
  dpd_contract244(&R1, &Z, &Z2, 0, 2, 1, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&Z);
  dpd_buf4_sort_axpy(&Z2, EOM_TMP0, pqsr, 10, 10, "Z(IA,jb)", -1.0);
  dpd_buf4_close(&Z2);

  psio_close(EOM_TMP1,0);
  psio_open(EOM_TMP1,PSIO_OPEN_NEW);

  /* - P(ij) P(ab) (Tjmbe Lnmfe) Rif Tna, term 29 */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 0, 10, 0, 0, "Z(JB,NI)");
  dpd_buf4_init(&I, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "VIAJB");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract424(&I, &R1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&I);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 11, 10, 11, 10, 0, "Z(AI,JB)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract244(&T1, &Z, &Z2, 0, 2, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_sort_axpy(&Z2, EOM_TMP0, qprs, 10, 10, "Z(IA,JB)", -1.0);
  dpd_buf4_close(&Z2);
  
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 0, 10, 0, 0, "Z(jb,ni)");
  dpd_buf4_init(&I, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "Viajb");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
  dpd_contract424(&I, &R1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&I);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 11, 10, 11, 10, 0, "Z(ai,jb)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract244(&T1, &Z, &Z2, 0, 2, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_sort_axpy(&Z2, EOM_TMP0, qprs, 10, 10, "Z(ia,jb)", -1.0);
  dpd_buf4_close(&Z2);

  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 0, 10, 0, 0, "Z(jb,NI)");
  dpd_buf4_init(&I, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "ViaJB");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract424(&I, &R1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&I);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 11, 10, 11, 10, 0, "Z(AI,jb)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract244(&T1, &Z, &Z2, 0, 2, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_sort_axpy(&Z2, EOM_TMP0, qprs, 10, 10, "Z(IA,jb)", -1.0);
  dpd_buf4_close(&Z2);

  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 0, 10, 0, 0, "Z(Ib,Nj)");
  dpd_buf4_init(&I, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "VIaJb");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
  dpd_contract424(&I, &R1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&I);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(Ib,Aj)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract244(&T1, &Z, &Z2, 0, 2, 1, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_sort_axpy(&Z2, EOM_TMP0, prsq, 10, 10, "Z(IA,jb)", +1.0);
  dpd_buf4_close(&Z2);

  psio_close(EOM_TMP1,0);
  psio_open(EOM_TMP1,PSIO_OPEN_NEW);

  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 0, 10, 0, 0, "Z(jA,nI)");
  dpd_buf4_init(&I, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "ViAjB");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract424(&I, &R1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&I);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 11, 10, 11, 10, 0, "Z(bI,jA)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract244(&T1, &Z, &Z2, 0, 2, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_sort(&Z2, EOM_TMP1, sqrp, 11, 10, "Z(AI,jb)");
  dpd_buf4_close(&Z2);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 11, 10, 11, 10, 0, "Z(AI,jb)");
  dpd_buf4_sort_axpy(&Z2, EOM_TMP0, qprs, 10, 10, "Z(IA,jb)", +1.0);
  dpd_buf4_close(&Z2);

  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 0, 10, 0, 0, "Z(IA,nj)");
  dpd_buf4_init(&I, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "VIAjb");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
  dpd_contract424(&I, &R1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&I);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(IA,bj)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract244(&T1, &Z, &Z2, 0, 2, 1, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_sort_axpy(&Z2, EOM_TMP0, pqsr, 10, 10, "Z(IA,jb)", -1.0);
  dpd_buf4_close(&Z2);

  psio_close(EOM_TMP1,0);
  psio_open(EOM_TMP1,PSIO_OPEN_NEW);

  /* - P(ij) P(ab) (Rjmbe Lnmfe) Tif Tna, term 21 */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 0, 10, 0, 0, "Z(JB,NI)");
  dpd_buf4_init(&I, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVOV");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&I, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&I);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 11, 10, 11, 10, 0, "Z(AI,JB)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract244(&T1, &Z, &Z2, 0, 2, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_sort_axpy(&Z2, EOM_TMP0, qprs, 10, 10, "Z(IA,JB)", -1.0);
  dpd_buf4_close(&Z2);
  
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 0, 10, 0, 0, "Z(jb,ni)");
  dpd_buf4_init(&I, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_ovov");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&I, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&I);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 11, 10, 11, 10, 0, "Z(ai,jb)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract244(&T1, &Z, &Z2, 0, 2, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_sort_axpy(&Z2, EOM_TMP0, qprs, 10, 10, "Z(ia,jb)", -1.0);
  dpd_buf4_close(&Z2);

  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 0, 10, 0, 0, "Z(jb,NI)");
  dpd_buf4_init(&I, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_ovOV");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&I, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&I);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 11, 10, 11, 10, 0, "Z(AI,jb)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract244(&T1, &Z, &Z2, 0, 2, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_sort_axpy(&Z2, EOM_TMP0, qprs, 10, 10, "Z(IA,jb)", -1.0);
  dpd_buf4_close(&Z2);

  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 0, 10, 0, 0, "Z(Ib,Nj)");
  dpd_buf4_init(&I, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OvOv");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&I, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&I);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(Ib,Aj)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract244(&T1, &Z, &Z2, 0, 2, 1, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_sort_axpy(&Z2, EOM_TMP0, prsq, 10, 10, "Z(IA,jb)", +1.0);
  dpd_buf4_close(&Z2);

  psio_close(EOM_TMP1,0);
  psio_open(EOM_TMP1,PSIO_OPEN_NEW);

  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 0, 10, 0, 0, "Z(jA,nI)");
  dpd_buf4_init(&I, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_oVoV");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&I, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&I);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 11, 10, 11, 10, 0, "Z(bI,jA)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract244(&T1, &Z, &Z2, 0, 2, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_sort(&Z2, EOM_TMP1, sqrp, 11, 10, "Z(AI,jb)");
  dpd_buf4_close(&Z2);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 11, 10, 11, 10, 0, "Z(AI,jb)");
  dpd_buf4_sort_axpy(&Z2, EOM_TMP0, qprs, 10, 10, "Z(IA,jb)", +1.0);
  dpd_buf4_close(&Z2);

  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 0, 10, 0, 0, "Z(IA,nj)");
  dpd_buf4_init(&I, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVov");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&I, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&I);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 10, 11, 10, 11, 0, "Z(IA,bj)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract244(&T1, &Z, &Z2, 0, 2, 1, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_sort_axpy(&Z2, EOM_TMP0, pqsr, 10, 10, "Z(IA,jb)", -1.0);
  dpd_buf4_close(&Z2);

  psio_close(EOM_TMP1,0);
  psio_open(EOM_TMP1,PSIO_OPEN_NEW);


  /* Now permute Z(IA,JB) and Z(ia,jb) and add them in along with Z(IA,jb) */
  dpd_buf4_init(&Z2, EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "Z(IA,JB)");
  dpd_buf4_sort(&Z2, EOM_TMP1, prqs, 0, 5, "Z(IJ,AB)");
  dpd_buf4_close(&Z2);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 5, 2, 7, 0, "GIJAB");
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z(IJ,AB)");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_sort(&Z, EOM_TMP1, qprs, 0, 5, "Z(JI,AB)");
  dpd_buf4_sort(&Z, EOM_TMP1, pqsr, 0, 5, "Z(IJ,BA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z(JI,AB)");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z(IJ,BA)");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_sort(&Z, EOM_TMP1, qprs, 0, 5, "Z(JI,BA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z(JI,BA)");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  dpd_buf4_init(&Z2, EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "Z(ia,jb)");
  dpd_buf4_sort(&Z2, EOM_TMP1, prqs, 0, 5, "Z(ij,ab)");
  dpd_buf4_close(&Z2);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 5, 2, 7, 0, "Gijab");
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z(ij,ab)");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_sort(&Z, EOM_TMP1, qprs, 0, 5, "Z(ji,ab)");
  dpd_buf4_sort(&Z, EOM_TMP1, pqsr, 0, 5, "Z(ij,ba)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z(ji,ab)");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z(ij,ba)");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_sort(&Z, EOM_TMP1, qprs, 0, 5, "Z(ji,ba)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z(ji,ba)");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 10, 10, 10, 10, 0, "Z(IA,jb)");
  dpd_buf4_sort_axpy(&Z, EOM_TMP0, prqs, 0, 5, "GIjAb", 1.0);
  dpd_buf4_close(&Z);

  /* add to ground state parts */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 2, 7, 2, 7, 0, "GIJAB");
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 2, 7, 2, 7, 0, "GIJAB");
  dpd_buf4_axpy(&G, &V, 0.5);
  dpd_buf4_close(&V);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 2, 7, 2, 7, 0, "Gijab");
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 2, 7, 2, 7, 0, "Gijab");
  dpd_buf4_axpy(&G, &V, 0.5);
  dpd_buf4_close(&V);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 5, 0, 5, 0, "GIjAb");
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 0, 5, 0, 5, 0, "GIjAb");
  dpd_buf4_axpy(&G, &V, 0.5);
  dpd_buf4_close(&V);
  dpd_buf4_close(&G);

  /*
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 2, 7, 2, 7, 0, "GIJAB");
  tval = dpd_buf4_dot_self(&V);
  dpd_buf4_close(&V);
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 2, 7, 2, 7, 0, "Gijab");
  tval += dpd_buf4_dot_self(&V);
  dpd_buf4_close(&V);
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 0, 5, 0, 5, 0, "GIjAb");
  tval += dpd_buf4_dot_self(&V);
  dpd_buf4_close(&V);
  fprintf(outfile,"<Gijab|Gijab> = %15.10lf\n", tval);
  */


  return;
}
