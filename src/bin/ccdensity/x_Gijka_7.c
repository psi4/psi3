#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/* This function computes 
   Gijka -= P(ij) LT_OO(k,i) * R(j,a) */

void x_Gijka_7(void) { 
  int h, nirreps, i, j, k, a, I, J, K, A, Isym, Jsym, Ksym, Asym, row, col;
  int II,JJ,IIsym,JJsym;
  int L_irr, R_irr, G_irr;
  dpdfile2 R1A, R1B, LTA, LTB;
  dpdbuf4 G;
 
  L_irr = params.L_irr; R_irr = params.R_irr; G_irr = params.G_irr;
  nirreps = moinfo.nirreps;

  /* open one-electron files for the nasty terms */
  dpd_file2_init(&R1A, CC_GR, R_irr, 0, 1, "RIA");
  dpd_file2_init(&R1B, CC_GR, R_irr, 0, 1, "Ria");
  dpd_file2_init(&LTA, EOM_TMP, L_irr, 0, 0, "LT_OO");
  dpd_file2_init(&LTB, EOM_TMP, L_irr, 0, 0, "LT_oo");
  dpd_file2_mat_init(&R1A);
  dpd_file2_mat_init(&R1B);
  dpd_file2_mat_init(&LTA);
  dpd_file2_mat_init(&LTB);
  dpd_file2_mat_rd(&R1A);
  dpd_file2_mat_rd(&R1B);
  dpd_file2_mat_rd(&LTA);
  dpd_file2_mat_rd(&LTB);

  /* rho_IJKA += - LT_OO(K,I)   R(J,A) + LT_OO(K,J)   R(I,A) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 10, 2, 10, 0, "GIJKA");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      j = G.params->roworb[h][row][1];
      I = LTA.params->colidx[i];  Isym = LTA.params->qsym[i];
      J = R1A.params->rowidx[j];  Jsym = R1A.params->psym[j];
      II = R1A.params->rowidx[i]; IIsym = R1A.params->psym[i];
      JJ = LTA.params->colidx[j]; JJsym = LTA.params->qsym[j];
      for(col=0; col < G.params->coltot[h]; col++) {
        k = G.params->colorb[h][col][0];
        a = G.params->colorb[h][col][1];
        K = LTA.params->rowidx[k]; Ksym = LTA.params->psym[k];
        A = R1A.params->colidx[a]; Asym = R1A.params->qsym[a];
        if( ((Ksym^Isym)==G_irr) && (Jsym==Asym))
          G.matrix[h][row][col] -=
            LTA.matrix[Ksym][K][I] * R1A.matrix[Jsym][J][A];
        if( ((Ksym^JJsym)==G_irr) && (IIsym==Asym))
          G.matrix[h][row][col] +=
            LTA.matrix[Ksym][K][JJ] * R1A.matrix[IIsym][II][A];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);

  /* rho_ijka += - LT_oo(k,i)   R(j,a) + LT_oo(k,j)   R(i,a) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 10, 2, 10, 0, "Gijka");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      j = G.params->roworb[h][row][1];
      I = LTB.params->colidx[i];  Isym = LTB.params->qsym[i];
      J = R1B.params->rowidx[j];  Jsym = R1B.params->psym[j];
      II = R1B.params->rowidx[i]; IIsym = R1B.params->psym[i];
      JJ = LTB.params->colidx[j]; JJsym = LTB.params->qsym[j];
      for(col=0; col < G.params->coltot[h]; col++) {
        k = G.params->colorb[h][col][0];
        a = G.params->colorb[h][col][1];
        K = LTB.params->rowidx[k]; Ksym = LTB.params->psym[k];
        A = R1B.params->colidx[a]; Asym = R1B.params->qsym[a];
        if( ((Ksym^Isym)==G_irr) && (Jsym==Asym))
          G.matrix[h][row][col] -=
            LTB.matrix[Ksym][K][I] * R1B.matrix[Jsym][J][A];
        if( ((Ksym^JJsym)==G_irr) && (IIsym==Asym))
          G.matrix[h][row][col] +=
            LTB.matrix[Ksym][K][JJ] * R1B.matrix[IIsym][II][A];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);

  /* rho_IjKa += - LT_OO(K,I)   R(j,a) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GIjKa");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      j = G.params->roworb[h][row][1];
      I = LTA.params->colidx[i];  Isym = LTA.params->qsym[i];
      J = R1B.params->rowidx[j];  Jsym = R1B.params->psym[j];
      for(col=0; col < G.params->coltot[h]; col++) {
        k = G.params->colorb[h][col][0];
        a = G.params->colorb[h][col][1];
        K = LTA.params->rowidx[k]; Ksym = LTA.params->psym[k];
        A = R1B.params->colidx[a]; Asym = R1B.params->qsym[a];
        if( ((Ksym^Isym)==G_irr) && (Jsym==Asym))
          G.matrix[h][row][col] -=
            LTA.matrix[Ksym][K][I] * R1B.matrix[Jsym][J][A];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);

  /* rho_iJkA += - LT_oo(k,i)   R(J,A) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GiJkA");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      j = G.params->roworb[h][row][1];
      I = LTB.params->colidx[i];  Isym = LTB.params->qsym[i];
      J = R1A.params->rowidx[j];  Jsym = R1A.params->psym[j];
      for(col=0; col < G.params->coltot[h]; col++) {
        k = G.params->colorb[h][col][0];
        a = G.params->colorb[h][col][1];
        K = LTB.params->rowidx[k]; Ksym = LTB.params->psym[k];
        A = R1A.params->colidx[a]; Asym = R1A.params->qsym[a];
        if( ((Ksym^Isym)==G_irr) && (Jsym==Asym))
          G.matrix[h][row][col] -=
            LTB.matrix[Ksym][K][I] * R1A.matrix[Jsym][J][A];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);

  dpd_file2_mat_close(&LTA);
  dpd_file2_mat_close(&LTB);
  dpd_file2_mat_close(&R1A);
  dpd_file2_mat_close(&R1B);

  dpd_file2_close(&LTA);
  dpd_file2_close(&LTB);
  dpd_file2_close(&R1A);
  dpd_file2_close(&R1B);

  return;
}
