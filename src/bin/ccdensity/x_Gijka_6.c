#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/* This function computes term 6,
   rho_ijka -= P(ij) lke rie tja or
   rho_ijka -= P(ij) LR1_OO(k,i) T(j,a)
 */

void x_Gijka_6(void) { 
  int h, nirreps, i, j, k, a, I, J, K, A, Isym, Jsym, Ksym, Asym, row, col;
  int II,JJ,IIsym,JJsym;
  int L_irr, R_irr, G_irr;
  dpdfile2 LR1A, LR1B, T1A, T1B;
  dpdbuf4 G;
 
  L_irr = params.L_irr; R_irr = params.R_irr; G_irr = params.G_irr;
  nirreps = moinfo.nirreps;

  /* open one-electron files for the nasty terms */
  dpd_file2_init(&LR1A, EOM_TMP, G_irr, 0, 0, "LR_OO");
  dpd_file2_init(&LR1B, EOM_TMP, G_irr, 0, 0, "LR_oo");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_init(&T1B, CC_OEI, 0, 0, 1, "tia");
  dpd_file2_mat_init(&T1A);   dpd_file2_mat_init(&T1B);
  dpd_file2_mat_init(&LR1A);   dpd_file2_mat_init(&LR1B);
  dpd_file2_mat_rd(&T1A);     dpd_file2_mat_rd(&T1B);
  dpd_file2_mat_rd(&LR1A);     dpd_file2_mat_rd(&LR1B);

  /* rho_IJKA += - LR1_OO(K,I) T(J,A) + LR1_OO(K,J) T(I,A) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 10, 2, 10, 0, "GIJKA");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      j = G.params->roworb[h][row][1];
      I = LR1A.params->colidx[i]; Isym = LR1A.params->qsym[i];
      J = T1A.params->rowidx[j]; Jsym = T1A.params->psym[j];
      II = T1A.params->rowidx[i]; IIsym = T1A.params->psym[i];
      JJ = LR1A.params->colidx[j]; JJsym = LR1A.params->qsym[j];
      for(col=0; col < G.params->coltot[h^G_irr]; col++) {
        k = G.params->colorb[h^G_irr][col][0];
        a = G.params->colorb[h^G_irr][col][1];
        K = LR1A.params->rowidx[k]; Ksym = LR1A.params->psym[k];
        A = T1A.params->colidx[a]; Asym = T1A.params->qsym[a];
        if( ((Ksym^Isym)==G_irr) && (Jsym==Asym))
          G.matrix[h][row][col] -=
            LR1A.matrix[Ksym][K][I] * T1A.matrix[Jsym][J][A];
        if( ((Ksym^JJsym)==G_irr) && (IIsym==Asym))
          G.matrix[h][row][col] +=
            LR1A.matrix[Ksym][K][JJ] * T1A.matrix[IIsym][II][A];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);

  /* rho_ijka += - LR1_oo(k,i) T(j,a) + LR1_oo(k,j) T(i,a) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 10, 2, 10, 0, "Gijka");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      j = G.params->roworb[h][row][1];
      I = LR1B.params->colidx[i]; Isym = LR1B.params->qsym[i];
      J = T1B.params->rowidx[j]; Jsym = T1B.params->psym[j];
      II = T1B.params->rowidx[i]; IIsym = T1B.params->psym[i];
      JJ = LR1B.params->colidx[j]; JJsym = LR1B.params->qsym[j];
      for(col=0; col < G.params->coltot[h^G_irr]; col++) {
        k = G.params->colorb[h^G_irr][col][0];
        a = G.params->colorb[h^G_irr][col][1];
        K = LR1B.params->rowidx[k]; Ksym = LR1B.params->psym[k];
        A = T1B.params->colidx[a]; Asym = T1B.params->qsym[a];
        if( ((Ksym^Isym)==G_irr) && (Jsym==Asym))
          G.matrix[h][row][col] -=
            LR1B.matrix[Ksym][K][I] * T1B.matrix[Jsym][J][A];
        if( ((Ksym^JJsym)==G_irr) && (IIsym==Asym))
          G.matrix[h][row][col] +=
            LR1B.matrix[Ksym][K][JJ] * T1B.matrix[IIsym][II][A];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);

  /* rho_IjKa += - LR1_OO(K,I) T(j,a) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GIjKa");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      j = G.params->roworb[h][row][1];
      I = LR1A.params->colidx[i]; Isym = LR1A.params->qsym[i];
      J = T1B.params->rowidx[j]; Jsym = T1B.params->psym[j];
      for(col=0; col < G.params->coltot[h^G_irr]; col++) {
        k = G.params->colorb[h^G_irr][col][0];
        a = G.params->colorb[h^G_irr][col][1];
        K = LR1A.params->rowidx[k]; Ksym = LR1A.params->psym[k];
        A = T1B.params->colidx[a]; Asym = T1B.params->qsym[a];
        if( ((Ksym^Isym)==G_irr) && (Jsym==Asym))
          G.matrix[h][row][col] -=
            LR1A.matrix[Ksym][K][I] * T1B.matrix[Jsym][J][A];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);

  /* rho_iJkA += - LR1_oo(k,i) T(J,A) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GiJkA");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      j = G.params->roworb[h][row][1];
      I = LR1B.params->colidx[i]; Isym = LR1B.params->qsym[i];
      J = T1A.params->rowidx[j]; Jsym = T1A.params->psym[j];
      for(col=0; col < G.params->coltot[h^G_irr]; col++) {
        k = G.params->colorb[h^G_irr][col][0];
        a = G.params->colorb[h^G_irr][col][1];
        K = LR1B.params->rowidx[k]; Ksym = LR1B.params->psym[k];
        A = T1A.params->colidx[a]; Asym = T1A.params->qsym[a];
        if( ((Ksym^Isym)==G_irr) && (Jsym==Asym))
          G.matrix[h][row][col] -=
            LR1B.matrix[Ksym][K][I] * T1A.matrix[Jsym][J][A];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);
  dpd_file2_mat_close(&LR1A);
  dpd_file2_mat_close(&LR1B);
  dpd_file2_close(&LR1A);
  dpd_file2_close(&LR1B);

  dpd_file2_mat_close(&T1A);
  dpd_file2_mat_close(&T1B);
  dpd_file2_close(&T1A);
  dpd_file2_close(&T1B);

  return;
}
