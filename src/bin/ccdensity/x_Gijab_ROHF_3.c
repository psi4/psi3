#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/* This function computes the following EOM Gijab terms
   P(ij) P(ab) [ Z(i,a) * R(j,b) ]
   Z(i,a) += Timae Lme, term 18
   Z(i,a) -= 0.5 (Lmnef Tmnea) Tif, term 32
   Z(i,a) -= 0.5 (Lmnef Tmief) Tna; term 34
   Z(i,a) += Lme Tma Tie ; term 19
 */

void x_Gijab_ROHF_3(void) { 
  int h, nirreps, row, col;
  int i,j,a,b;
  int I1, I2, I3, I4, J1, J2, J3, J4, A1, A2, A3, A4, B1, B2, B3, B4;
  int I1sym, I2sym, I3sym, I4sym, J1sym, J2sym, J3sym, J4sym;
  int A1sym, A2sym, A3sym, A4sym, B1sym, B2sym, B3sym, B4sym;
  int L_irr, R_irr, G_irr;
  dpdfile2 L1T2A, L1T2B, T1A, T1B, Z1A, Z1B, I1A, I1B, R1A, R1B;
  dpdbuf4 G, I, Z, Z2, T, L;
 
  L_irr = params.L_irr; R_irr = params.R_irr; G_irr = params.G_irr;
  nirreps = moinfo.nirreps;

  /* Z(I,A) += L1T2_OV, term 18 */
  dpd_file2_init(&L1T2A, EOM_TMP, L_irr, 0, 1, "L1T2_OV");
  dpd_file2_init(&L1T2B, EOM_TMP, L_irr, 0, 1, "L1T2_ov");
  dpd_file2_copy(&L1T2A, EOM_TMP1, "ZIA");
  dpd_file2_copy(&L1T2B, EOM_TMP1, "Zia");
  dpd_file2_close(&L1T2A);
  dpd_file2_close(&L1T2B);

  dpd_file2_init(&Z1A, EOM_TMP1, L_irr, 0, 1, "ZIA");
  dpd_file2_init(&Z1B, EOM_TMP1, L_irr, 0, 1, "Zia");

  /* Z(I,A) -= 0.5 Lmnef Tmnea Tif, term 32 */
  dpd_file2_init(&I1A, EOM_TMP, L_irr, 1, 1, "LT2_VV");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract222(&T1A, &I1A, &Z1A, 0, 1, -1.0, 1.0);
  dpd_file2_close(&T1A);
  dpd_file2_close(&I1A);
  dpd_file2_init(&I1B, EOM_TMP, L_irr, 1, 1, "LT2_vv");
  dpd_file2_init(&T1B, CC_OEI, 0, 0, 1, "tia");
  dpd_contract222(&T1B, &I1B, &Z1B, 0, 1, -1.0, 1.0);
  dpd_file2_close(&T1B);
  dpd_file2_close(&I1B);

  /* Z(i,a) -= 0.5 (Lmnef Tmief ) Tna; term 34 */
  dpd_file2_init(&I1A, EOM_TMP, L_irr, 0, 0, "LT2_OO");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract222(&I1A, &T1A, &Z1A, 1, 1, -1.0, 1.0);
  dpd_file2_close(&T1A);
  dpd_file2_close(&I1A);
  dpd_file2_init(&I1B, EOM_TMP, L_irr, 0, 0, "LT2_oo");
  dpd_file2_init(&T1B, CC_OEI, 0, 0, 1, "tia");
  dpd_contract222(&I1B, &T1B, &Z1B, 1, 1, -1.0, 1.0);
  dpd_file2_close(&T1B);
  dpd_file2_close(&I1B);

  /* Z(i,a) += Lme Tma Tie ; term 19  */
  dpd_file2_init(&I1A, EOM_TMP, L_irr, 1, 1, "LT1_VV");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract222(&T1A, &I1A, &Z1A, 0, 1, -1.0, 1.0);
  dpd_file2_close(&T1A);
  dpd_file2_close(&I1A);
  dpd_file2_init(&I1B, EOM_TMP, L_irr, 1, 1, "LT1_vv");
  dpd_file2_init(&T1B, CC_OEI, 0, 0, 1, "tia");
  dpd_contract222(&T1B, &I1B, &Z1B, 0, 1, -1.0, 1.0);
  dpd_file2_close(&T1B);
  dpd_file2_close(&I1B);

  /* open one-electron files for the nasty terms */
  dpd_file2_mat_init(&Z1A);   dpd_file2_mat_init(&Z1B);
  dpd_file2_mat_rd(&Z1A);     dpd_file2_mat_rd(&Z1B);

  dpd_file2_init(&R1A, CC_GR, R_irr, 0, 1, "RIA");
  dpd_file2_init(&R1B, CC_GR, R_irr, 0, 1, "Ria");
  dpd_file2_mat_init(&R1A);   dpd_file2_mat_init(&R1B);
  dpd_file2_mat_rd(&R1A);     dpd_file2_mat_rd(&R1B);

  /* + Z(I,A) R(J,B) */
  /* - Z(I,B) R(J,A) */
  /* + R(I,A) Z(J,B) */
  /* - R(I,B) Z(J,A) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 2, 7, 2, 7, 0, "GIJAB");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      j = G.params->roworb[h][row][1];
      I1 = Z1A.params->rowidx[i]; I1sym = Z1A.params->psym[i];
      I2 = I1; I2sym = I1sym;
      I3 = R1A.params->rowidx[i]; I3sym = R1A.params->psym[i];
      I4 = I3; I4sym = I3sym;
      J1 = R1A.params->rowidx[j]; J1sym = R1A.params->psym[j];
      J2 = J1; J2sym=J1sym;
      J3 = Z1A.params->rowidx[j]; J3sym = Z1A.params->psym[j];
      J4 = J3; J4sym=J3sym;
      for(col=0; col < G.params->coltot[h]; col++) {
        a = G.params->colorb[h][col][0];
        b = G.params->colorb[h][col][1];
        A1 = Z1A.params->colidx[a]; A1sym = Z1A.params->qsym[a];
        A4 = A1; A4sym = A1sym;
        A2 = R1A.params->colidx[a]; A2sym = R1A.params->qsym[a];
        A3 = A2; A3sym = A2sym;
        B1 = R1A.params->colidx[b]; B1sym = R1A.params->qsym[b];
        B4 = B1; B4sym = B1sym;
        B2 = Z1A.params->colidx[b]; B2sym = Z1A.params->qsym[b];
        B3 = B2; B3sym = B2sym;
        /* + Z(I,A) R(J,B) */
        if ( ((I1sym^A1sym)==L_irr) && ((J1sym^B1sym)==R_irr) )
          G.matrix[h][row][col] +=
            Z1A.matrix[I1sym][I1][A1] * R1A.matrix[J1sym][J1][B1];
        /* - Z(I,B) R(J,A) */
        if ( ((I2sym^B2sym)==L_irr) && ((J2sym^A2sym)==R_irr) )
          G.matrix[h][row][col] -=
            Z1A.matrix[I2sym][I2][B2] * R1A.matrix[J2sym][J2][A2];
        /* + R(I,A) Z(J,B) */
        if ( ((J3sym^B3sym)==L_irr) && ((I3sym^A3sym)==R_irr) )
          G.matrix[h][row][col] +=
            Z1A.matrix[J3sym][J3][B3] * R1A.matrix[I3sym][I3][A3];
        /* - R(I,B) Z(J,A) */
        if ( ((J4sym^A4sym)==L_irr) && ((I4sym^B4sym)==R_irr) )
          G.matrix[h][row][col] -=
            Z1A.matrix[J4sym][J4][A4] * R1A.matrix[I4sym][I4][B4];
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
      I3 = R1B.params->rowidx[i]; I3sym = R1B.params->psym[i];
      I4 = I3; I4sym = I3sym;
      J1 = R1B.params->rowidx[j]; J1sym = R1B.params->psym[j];
      J2 = J1; J2sym=J1sym;
      J3 = Z1B.params->rowidx[j]; J3sym = Z1B.params->psym[j];
      J4 = J3; J4sym=J3sym;
      for(col=0; col < G.params->coltot[h]; col++) {
        a = G.params->colorb[h][col][0];
        b = G.params->colorb[h][col][1];
        A1 = Z1B.params->colidx[a]; A1sym = Z1B.params->qsym[a];
        A4 = A1; A4sym = A1sym;
        A2 = R1B.params->colidx[a]; A2sym = R1B.params->qsym[a];
        A3 = A2; A3sym = A2sym;
        B1 = R1B.params->colidx[b]; B1sym = R1B.params->qsym[b];
        B4 = B1; B4sym = B1sym;
        B2 = Z1B.params->colidx[b]; B2sym = Z1B.params->qsym[b];
        B3 = B2; B3sym = B2sym;
        /* + Z(i,a) R(j,b) */
        if ( ((I1sym^A1sym)==L_irr) && ((J1sym^B1sym)==R_irr ) )
          G.matrix[h][row][col] +=
            Z1B.matrix[I1sym][I1][A1] * R1B.matrix[J1sym][J1][B1];
        /* - Z(i,b) R(j,a) */
        if ( ((I2sym^B2sym)==L_irr) && ((J2sym^A2sym)==R_irr ) )
          G.matrix[h][row][col] -=
            Z1B.matrix[I2sym][I2][B2] * R1B.matrix[J2sym][J2][A2];
        /* + R(i,a) Z(j,b) */
        if ( ((J3sym^B3sym)==L_irr) && ((I3sym^A3sym)==R_irr ) )
          G.matrix[h][row][col] +=
            Z1B.matrix[J3sym][J3][B3] * R1B.matrix[I3sym][I3][A3];
        /* - R(i,b) Z(j,a) */
        if ( ((J4sym^A4sym)==L_irr) && ((I4sym^B4sym)==R_irr ) )
          G.matrix[h][row][col] -=
            Z1B.matrix[J4sym][J4][A4] * R1B.matrix[I4sym][I4][B4];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);

  /* + Z(I,A) R(j,b) */
  /* + R(I,A) Z(j,b) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 5, 0, 5, 0, "GIjAb");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      j = G.params->roworb[h][row][1];
      I1 = Z1A.params->rowidx[i]; I1sym = Z1A.params->psym[i];
      I2 = R1A.params->rowidx[i]; I2sym = R1A.params->psym[i];
      J1 = R1B.params->rowidx[j]; J1sym = R1B.params->psym[j];
      J2 = Z1B.params->rowidx[j]; J2sym = Z1B.params->psym[j];
      for(col=0; col < G.params->coltot[h]; col++) {
        a = G.params->colorb[h][col][0];
        b = G.params->colorb[h][col][1];
        A1 = Z1A.params->colidx[a]; A1sym = Z1A.params->qsym[a];
        A2 = R1A.params->colidx[a]; A2sym = R1A.params->qsym[a];
        B1 = R1B.params->colidx[b]; B1sym = R1B.params->qsym[b];
        B2 = Z1B.params->colidx[b]; B2sym = Z1B.params->qsym[b];
        /* + Z(I,A) R(j,b) */
        if ( ((I1sym^A1sym)==L_irr) && ((J1sym^B1sym)==R_irr ) )
          G.matrix[h][row][col] +=
            Z1A.matrix[I1sym][I1][A1] * R1B.matrix[J1sym][J1][B1];
        /* + R(I,A) Z(j,b) */
        if ( ((J2sym^B2sym)==L_irr) && ((I2sym^A2sym)==R_irr ) )
          G.matrix[h][row][col] +=
            R1A.matrix[I2sym][I2][A2] * Z1B.matrix[J2sym][J2][B2];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);

  dpd_file2_mat_close(&R1A); dpd_file2_mat_close(&R1B);
  dpd_file2_close(&R1A); dpd_file2_close(&R1B);

  dpd_file2_mat_close(&Z1A); dpd_file2_mat_close(&Z1B);
  dpd_file2_close(&Z1A); dpd_file2_close(&Z1B);

  psio_close(EOM_TMP1,0);
  psio_open(EOM_TMP1,PSIO_OPEN_NEW);
  return;
}
