#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/* This function computes the following EOM Gijab terms
   Z(i,a) += Rimae Lme ; term 3
   Z(i,a) -= 0.5 Lmnef Tmnea Rif ; term 22
   Z(i,a) -= 0.5 Lmnef Timef Rna; term 23
   Z(i,a) += lmnef rme tinaf; term 33 
   P(ij) P(ab) [ Z(i,a) * T(j,b) ]
 */

void x_Gijab_ROHF_2(void) { 
  int h, nirreps, row, col;
  int i,j,a,b;
  int I1, I2, I3, I4, J1, J2, J3, J4, A1, A2, A3, A4, B1, B2, B3, B4;
  int I1sym, I2sym, I3sym, I4sym, J1sym, J2sym, J3sym, J4sym;
  int A1sym, A2sym, A3sym, A4sym, B1sym, B2sym, B3sym, B4sym;
  int L_irr, R_irr, G_irr;
  dpdfile2 L1R2A, L1R2B, T1A, T1B, Z1A, Z1B, I1A, I1B, R1A, R1B;
  dpdbuf4 G, I, Z, Z2, T, L;
 
  L_irr = params.L_irr; R_irr = params.R_irr; G_irr = params.G_irr;
  nirreps = moinfo.nirreps;

  /* Z(I,A) += L1R2_OV, term 3 */
  dpd_file2_init(&L1R2A, EOM_TMP, G_irr, 0, 1, "L1R2_OV");
  dpd_file2_init(&L1R2B, EOM_TMP, G_irr, 0, 1, "L1R2_ov");
  dpd_file2_copy(&L1R2A, EOM_TMP1, "ZIA");
  dpd_file2_copy(&L1R2B, EOM_TMP1, "Zia");
  dpd_file2_close(&L1R2A);
  dpd_file2_close(&L1R2B);

  dpd_file2_init(&Z1A, EOM_TMP1, G_irr, 0, 1, "ZIA");
  dpd_file2_init(&Z1B, EOM_TMP1, G_irr, 0, 1, "Zia");

  /* Z(I,A) += 0.5 (lmnef tmnea) rif, term 22 */
  dpd_file2_init(&I1A, EOM_TMP, L_irr, 1, 1, "LT2_VV");
  dpd_file2_init(&R1A, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract222(&R1A, &I1A, &Z1A, 0, 1, -1.0, 1.0);
  dpd_file2_close(&R1A);
  dpd_file2_close(&I1A);
  dpd_file2_init(&I1B, EOM_TMP, L_irr, 1, 1, "LT2_vv");
  dpd_file2_init(&R1B, CC_GR, R_irr, 0, 1, "Ria");
  dpd_contract222(&R1B, &I1B, &Z1B, 0, 1, -1.0, 1.0);
  dpd_file2_close(&R1B);
  dpd_file2_close(&I1B);

  /* Z(i,a) -= 0.5 (timef lnmef) rna; term 23 */
  dpd_file2_init(&I1A, EOM_TMP, L_irr, 0, 0, "LT2_OO");
  dpd_file2_init(&R1A, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract222(&I1A, &R1A, &Z1A, 1, 1, -1.0, 1.0);
  dpd_file2_close(&R1A);
  dpd_file2_close(&I1A);
  dpd_file2_init(&I1B, EOM_TMP, L_irr, 0, 0, "LT2_oo");
  dpd_file2_init(&R1B, CC_GR, R_irr, 0, 1, "Ria");
  dpd_contract222(&I1B, &R1B, &Z1B, 1, 1, -1.0, 1.0);
  dpd_file2_close(&R1B);
  dpd_file2_close(&I1B);

  /* Z(i,a) += lmnef rme tinaf; term 33  */
  dpd_file2_init(&I1A, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
  dpd_file2_init(&I1B, EOM_TMP, G_irr, 0, 1, "L2R1_ov");
  dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
  dpd_dot24(&I1A, &T, &Z1A, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&T);
  dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_dot24(&I1B, &T, &Z1A, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&T);
  dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tijab");
  dpd_dot24(&I1B, &T, &Z1B, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&T);
  dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
  dpd_dot24(&I1A, &T, &Z1B, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&T);
  dpd_file2_close(&I1A);
  dpd_file2_close(&I1B);


  /* open one-electron files for the nasty permutations */
  dpd_file2_mat_init(&Z1A);   dpd_file2_mat_init(&Z1B);
  dpd_file2_mat_rd(&Z1A);     dpd_file2_mat_rd(&Z1B);

  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_init(&T1B, CC_OEI, 0, 0, 1, "tia");
  dpd_file2_mat_init(&T1A);   dpd_file2_mat_init(&T1B);
  dpd_file2_mat_rd(&T1A);     dpd_file2_mat_rd(&T1B);

  /* + Z(I,A) T(J,B) */
  /* - Z(I,B) T(J,A) */
  /* + T(I,A) Z(J,B) */
  /* - T(I,B) Z(J,A) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 2, 7, 2, 7, 0, "GIJAB");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      j = G.params->roworb[h][row][1];
      I1 = Z1A.params->rowidx[i]; I1sym = Z1A.params->psym[i];
      I2 = I1; I2sym = I1sym;
      I3 = T1A.params->rowidx[i]; I3sym = T1A.params->psym[i];
      I4 = I3; I4sym = I3sym;
      J1 = T1A.params->rowidx[j]; J1sym = T1A.params->psym[j];
      J2 = J1; J2sym=J1sym;
      J3 = Z1A.params->rowidx[j]; J3sym = Z1A.params->psym[j];
      J4 = J3; J4sym=J3sym;
      for(col=0; col < G.params->coltot[h]; col++) {
        a = G.params->colorb[h][col][0];
        b = G.params->colorb[h][col][1];
        A1 = Z1A.params->colidx[a]; A1sym = Z1A.params->qsym[a];
        A4 = A1; A4sym = A1sym;
        A2 = T1A.params->colidx[a]; A2sym = T1A.params->qsym[a];
        A3 = A2; A3sym = A2sym;
        B1 = T1A.params->colidx[b]; B1sym = T1A.params->qsym[b];
        B4 = B1; B4sym = B1sym;
        B2 = Z1A.params->colidx[b]; B2sym = Z1A.params->qsym[b];
        B3 = B2; B3sym = B2sym;
        /* + Z(I,A) T(J,B) */
        if ( ((I1sym^A1sym)==G_irr) && (J1sym==B1sym) )
          G.matrix[h][row][col] +=
            Z1A.matrix[I1sym][I1][A1] * T1A.matrix[J1sym][J1][B1];
        /* - Z(I,B) T(J,A) */
        if ( ((I2sym^B2sym)==G_irr) && (J2sym==A2sym) )
          G.matrix[h][row][col] -=
            Z1A.matrix[I2sym][I2][B2] * T1A.matrix[J2sym][J2][A2];
        /* + T(I,A) Z(J,B) */
        if ( ((J3sym^B3sym)==G_irr) && (I3sym==A3sym) )
          G.matrix[h][row][col] +=
            Z1A.matrix[J3sym][J3][B3] * T1A.matrix[I3sym][I3][A3];
        /* - T(I,B) Z(J,A) */
        if ( ((J4sym^A4sym)==G_irr) && (I4sym==B4sym) )
          G.matrix[h][row][col] -=
            Z1A.matrix[J4sym][J4][A4] * T1A.matrix[I4sym][I4][B4];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, EOM_TMP0, G_irr, 2, 7, 2, 7, 0, "Gijab");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      j = G.params->roworb[h][row][1];
      I1 = Z1B.params->rowidx[i]; I1sym = Z1B.params->psym[i];
      I2 = I1; I2sym = I1sym;
      I3 = T1B.params->rowidx[i]; I3sym = T1B.params->psym[i];
      I4 = I3; I4sym = I3sym;
      J1 = T1B.params->rowidx[j]; J1sym = T1B.params->psym[j];
      J2 = J1; J2sym=J1sym;
      J3 = Z1B.params->rowidx[j]; J3sym = Z1B.params->psym[j];
      J4 = J3; J4sym=J3sym;
      for(col=0; col < G.params->coltot[h]; col++) {
        a = G.params->colorb[h][col][0];
        b = G.params->colorb[h][col][1];
        A1 = Z1B.params->colidx[a]; A1sym = Z1B.params->qsym[a];
        A4 = A1; A4sym = A1sym;
        A2 = T1B.params->colidx[a]; A2sym = T1B.params->qsym[a];
        A3 = A2; A3sym = A2sym;
        B1 = T1B.params->colidx[b]; B1sym = T1B.params->qsym[b];
        B4 = B1; B4sym = B1sym;
        B2 = Z1B.params->colidx[b]; B2sym = Z1B.params->qsym[b];
        B3 = B2; B3sym = B2sym;
        /* + Z(i,a) T(j,b) */
        if ( ((I1sym^A1sym)==G_irr) && (J1sym==B1sym) )
          G.matrix[h][row][col] +=
            Z1B.matrix[I1sym][I1][A1] * T1B.matrix[J1sym][J1][B1];
        /* - Z(i,b) T(j,a) */
        if ( ((I2sym^B2sym)==G_irr) && (J2sym==A2sym) )
          G.matrix[h][row][col] -=
            Z1B.matrix[I2sym][I2][B2] * T1B.matrix[J2sym][J2][A2];
        /* + T(i,a) Z(j,b) */
        if ( ((J3sym^B3sym)==G_irr) && (I3sym==A3sym) )
          G.matrix[h][row][col] +=
            Z1B.matrix[J3sym][J3][B3] * T1B.matrix[I3sym][I3][A3];
        /* - T(i,b) Z(j,a) */
        if ( ((J4sym^A4sym)==G_irr) && (I4sym==B4sym) )
          G.matrix[h][row][col] -=
            Z1B.matrix[J4sym][J4][A4] * T1B.matrix[I4sym][I4][B4];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);

  /* + Z(I,A) T(j,b) */
  /* + T(I,A) Z(j,b) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 5, 0, 5, 0, "GIjAb");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      j = G.params->roworb[h][row][1];
      I1 = Z1A.params->rowidx[i]; I1sym = Z1A.params->psym[i];
      I2 = T1A.params->rowidx[i]; I2sym = T1A.params->psym[i];
      J1 = T1B.params->rowidx[j]; J1sym = T1B.params->psym[j];
      J2 = Z1B.params->rowidx[j]; J2sym = Z1B.params->psym[j];
      for(col=0; col < G.params->coltot[h]; col++) {
        a = G.params->colorb[h][col][0];
        b = G.params->colorb[h][col][1];
        A1 = Z1A.params->colidx[a]; A1sym = Z1A.params->qsym[a];
        A2 = T1A.params->colidx[a]; A2sym = T1A.params->qsym[a];
        B1 = T1B.params->colidx[b]; B1sym = T1B.params->qsym[b];
        B2 = Z1B.params->colidx[b]; B2sym = Z1B.params->qsym[b];
        /* + Z(I,A) T(j,b) */
        if ( ((I1sym^A1sym)==G_irr) && (J1sym==B1sym) )
          G.matrix[h][row][col] +=
            Z1A.matrix[I1sym][I1][A1] * T1B.matrix[J1sym][J1][B1];
        /* + T(I,A) Z(j,b) */
        if ( ((J2sym^B2sym)==G_irr) && (I2sym==A2sym) )
          G.matrix[h][row][col] +=
            T1A.matrix[I2sym][I2][A2] * Z1B.matrix[J2sym][J2][B2];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);

  dpd_file2_mat_close(&T1A); dpd_file2_mat_close(&T1B);
  dpd_file2_close(&T1A); dpd_file2_close(&T1B);

  dpd_file2_mat_close(&Z1A); dpd_file2_mat_close(&Z1B);
  dpd_file2_close(&Z1A); dpd_file2_close(&Z1B);

  psio_close(EOM_TMP1,0);
  psio_open(EOM_TMP1,PSIO_OPEN_NEW);
  return;
}
