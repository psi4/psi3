#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/* This function computes term 6,
   rho_ciab += P(ab) LR1_VV(c,a) T(i,b)
*/

void x_Gciab_6(void) { 
  int h, nirreps, c, i, a, b, C, I, A, B, Csym, Isym, Asym, Bsym, row, col;
  int AA, BB, AAsym, BBsym;
  int L_irr, R_irr, G_irr;
  dpdfile2 LR1A, LR1B, T1A, T1B;
  dpdbuf4 G;
 
  L_irr = params.L_irr; R_irr = params.R_irr; G_irr = params.G_irr;
  nirreps = moinfo.nirreps;

  /* open one-electron files for the nasty terms */
  dpd_file2_init(&LR1A, EOM_TMP, G_irr, 1, 1, "LR_VV");
  dpd_file2_init(&LR1B, EOM_TMP, G_irr, 1, 1, "LR_vv");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_init(&T1B, CC_OEI, 0, 0, 1, "tia");
  dpd_file2_mat_init(&T1A);   dpd_file2_mat_init(&T1B);
  dpd_file2_mat_init(&LR1A);   dpd_file2_mat_init(&LR1B);
  dpd_file2_mat_rd(&T1A);     dpd_file2_mat_rd(&T1B);
  dpd_file2_mat_rd(&LR1A);     dpd_file2_mat_rd(&LR1B);

  /* rho_CIAB += LR1_VV(C,A) T(I,B) - LR1_VV(C,B) T(I,A) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 7, 0, "GCIAB");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      c = G.params->roworb[h][row][0];
      i = G.params->roworb[h][row][1];
      C = LR1A.params->rowidx[c]; Csym = LR1A.params->psym[c];
      I = T1A.params->rowidx[i]; Isym = T1A.params->psym[i];
      for(col=0; col < G.params->coltot[h]; col++) {
        a = G.params->colorb[h^G_irr][col][0];
        b = G.params->colorb[h^G_irr][col][1];
        A = LR1A.params->colidx[a]; Asym = LR1A.params->qsym[a];
        B = T1A.params->colidx[b]; Bsym = T1A.params->qsym[b];
        AA = T1A.params->colidx[a]; AAsym = T1A.params->qsym[a];
        BB = LR1A.params->colidx[b]; BBsym = LR1A.params->qsym[b];
        if( ((Csym^Asym)==G_irr) && (Isym==Bsym))
          G.matrix[h][row][col] +=
            LR1A.matrix[Csym][C][A] * T1A.matrix[Isym][I][B];
        if( ((Csym^BBsym)==G_irr) && (Isym==AAsym))
          G.matrix[h][row][col] -=
            LR1A.matrix[Csym][C][BB] * T1A.matrix[Isym][I][AA];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);
  /* rho_ciab += LR1_vv(c,a) T(i,b) - LR1_vv(c,b) T(i,a) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 7, 0, "Gciab");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      c = G.params->roworb[h][row][0];
      i = G.params->roworb[h][row][1];
      C = LR1B.params->rowidx[c]; Csym = LR1B.params->psym[c];
      I = T1B.params->rowidx[i]; Isym = T1B.params->psym[i];
      for(col=0; col < G.params->coltot[h^G_irr]; col++) {
        a = G.params->colorb[h^G_irr][col][0];
        b = G.params->colorb[h^G_irr][col][1];
        A = LR1B.params->colidx[a]; Asym = LR1B.params->qsym[a];
        B = T1B.params->colidx[b]; Bsym = T1B.params->qsym[b];
        AA = T1B.params->colidx[a]; AAsym = T1B.params->qsym[a];
        BB = LR1B.params->colidx[b]; BBsym = LR1B.params->qsym[b];
        if( ((Csym^Asym)==G_irr) && (Isym==Bsym))
          G.matrix[h][row][col] +=
            LR1B.matrix[Csym][C][A] * T1B.matrix[Isym][I][B];
        if( ((Csym^BBsym)==G_irr) && (Isym==AAsym))
          G.matrix[h][row][col] -=
            LR1B.matrix[Csym][C][BB] * T1B.matrix[Isym][I][AA];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);
  /* rho_CiAb += LR1_VV(C,A) T(i,b) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GCiAb");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      c = G.params->roworb[h][row][0];
      i = G.params->roworb[h][row][1];
      C = LR1A.params->rowidx[c]; Csym = LR1A.params->psym[c];
      I = T1B.params->rowidx[i]; Isym = T1B.params->psym[i];
      for(col=0; col < G.params->coltot[h^G_irr]; col++) {
        a = G.params->colorb[h^G_irr][col][0];
        b = G.params->colorb[h^G_irr][col][1];
        A = LR1A.params->colidx[a]; Asym = LR1A.params->qsym[a];
        B = T1B.params->colidx[b]; Bsym = T1B.params->qsym[b];
        if( ((Csym^Asym)==G_irr) && (Isym==Bsym))
          G.matrix[h][row][col] +=
            LR1A.matrix[Csym][C][A] * T1B.matrix[Isym][I][B];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);
  /* rho_cIaB += LR1_vv(c,a) T(I,B) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GcIaB");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      c = G.params->roworb[h][row][0];
      i = G.params->roworb[h][row][1];
      C = LR1B.params->rowidx[c]; Csym = LR1B.params->psym[c];
      I = T1A.params->rowidx[i]; Isym = T1A.params->psym[i];
      for(col=0; col < G.params->coltot[h^G_irr]; col++) {
        a = G.params->colorb[h^G_irr][col][0];
        b = G.params->colorb[h^G_irr][col][1];
        A = LR1B.params->colidx[a]; Asym = LR1B.params->qsym[a];
        B = T1A.params->colidx[b]; Bsym = T1A.params->qsym[b];
        if( ((Csym^Asym)==G_irr) && (Isym==Bsym))
          G.matrix[h][row][col] +=
            LR1B.matrix[Csym][C][A] * T1A.matrix[Isym][I][B];
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
