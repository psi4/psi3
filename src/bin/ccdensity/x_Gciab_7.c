#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/* This function computes term 7,
   rho_ciab += P(ab) LT_VV(c,a) R(i,b)
*/

void x_Gciab_7(void) { 
  int h, nirreps, c, i, a, b, C, I, A, B, Csym, Isym, Asym, Bsym, row, col;
  int AA, BB, AAsym, BBsym;
  int L_irr, R_irr, G_irr;
  dpdfile2 LT1A, LT1B, R1A, R1B;
  dpdbuf4 G;
 
  L_irr = params.L_irr; R_irr = params.R_irr; G_irr = params.G_irr;
  nirreps = moinfo.nirreps;

  /* open one-electron files for the nasty terms */
  dpd_file2_init(&LT1A, EOM_TMP, L_irr, 1, 1, "LT_VV");
  dpd_file2_init(&LT1B, EOM_TMP, L_irr, 1, 1, "LT_vv");
  dpd_file2_init(&R1A, CC_GR, R_irr, 0, 1, "RIA");
  dpd_file2_init(&R1B, CC_GR, R_irr, 0, 1, "Ria");
  dpd_file2_mat_init(&R1A);   dpd_file2_mat_init(&R1B);
  dpd_file2_mat_init(&LT1A);   dpd_file2_mat_init(&LT1B);
  dpd_file2_mat_rd(&R1A);     dpd_file2_mat_rd(&R1B);
  dpd_file2_mat_rd(&LT1A);     dpd_file2_mat_rd(&LT1B);

  /* rho_CIAB += LT_VV(C,A) R(I,B) - LT_VV(C,B) R(I,A) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 7, 0, "GCIAB");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      c = G.params->roworb[h][row][0];
      i = G.params->roworb[h][row][1];
      C = LT1A.params->rowidx[c]; Csym = LT1A.params->psym[c];
      I = R1A.params->rowidx[i]; Isym = R1A.params->psym[i];
      for(col=0; col < G.params->coltot[h]; col++) {
        a = G.params->colorb[h][col][0];
        b = G.params->colorb[h][col][1];
        A = LT1A.params->colidx[a]; Asym = LT1A.params->qsym[a];
        B = R1A.params->colidx[b]; Bsym = R1A.params->qsym[b];
        AA = R1A.params->colidx[a]; AAsym = R1A.params->qsym[a];
        BB = LT1A.params->colidx[b]; BBsym = LT1A.params->qsym[b];
        if( ((Csym^Asym)==G_irr) && (Isym==Bsym))
          G.matrix[h][row][col] +=
            LT1A.matrix[Csym][C][A] * R1A.matrix[Isym][I][B];
        if( ((Csym^BBsym)==G_irr) && (Isym==AAsym))
          G.matrix[h][row][col] -=
            LT1A.matrix[Csym][C][BB] * R1A.matrix[Isym][I][AA];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);
  /* rho_ciab += LT_vv(c,a) R(i,b) - LT_vv(c,b) R(i,a) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 7, 0, "Gciab");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      c = G.params->roworb[h][row][0];
      i = G.params->roworb[h][row][1];
      C = LT1B.params->rowidx[c]; Csym = LT1B.params->psym[c];
      I = R1B.params->rowidx[i]; Isym = R1B.params->psym[i];
      for(col=0; col < G.params->coltot[h]; col++) {
        a = G.params->colorb[h][col][0];
        b = G.params->colorb[h][col][1];
        A = LT1B.params->colidx[a]; Asym = LT1B.params->qsym[a];
        B = R1B.params->colidx[b]; Bsym = R1B.params->qsym[b];
        AA = R1B.params->colidx[a]; AAsym = R1B.params->qsym[a];
        BB = LT1B.params->colidx[b]; BBsym = LT1B.params->qsym[b];
        if( ((Csym^Asym)==G_irr) && (Isym==Bsym))
          G.matrix[h][row][col] +=
            LT1B.matrix[Csym][C][A] * R1B.matrix[Isym][I][B];
        if( ((Csym^BBsym)==G_irr) && (Isym==AAsym))
          G.matrix[h][row][col] -=
            LT1B.matrix[Csym][C][BB] * R1B.matrix[Isym][I][AA];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);
  /* rho_CiAb += LT_VV(C,A) R(i,b) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GCiAb");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      c = G.params->roworb[h][row][0];
      i = G.params->roworb[h][row][1];
      C = LT1A.params->rowidx[c]; Csym = LT1A.params->psym[c];
      I = R1B.params->rowidx[i]; Isym = R1B.params->psym[i];
      for(col=0; col < G.params->coltot[h]; col++) {
        a = G.params->colorb[h][col][0];
        b = G.params->colorb[h][col][1];
        A = LT1A.params->colidx[a]; Asym = LT1A.params->qsym[a];
        B = R1B.params->colidx[b]; Bsym = R1B.params->qsym[b];
        if( ((Csym^Asym)==G_irr) && (Isym==Bsym))
          G.matrix[h][row][col] +=
            LT1A.matrix[Csym][C][A] * R1B.matrix[Isym][I][B];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);
  /* rho_cIaB += LT_vv(c,a) R(I,B) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GcIaB");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      c = G.params->roworb[h][row][0];
      i = G.params->roworb[h][row][1];
      C = LT1B.params->rowidx[c]; Csym = LT1B.params->psym[c];
      I = R1A.params->rowidx[i]; Isym = R1A.params->psym[i];
      for(col=0; col < G.params->coltot[h]; col++) {
        a = G.params->colorb[h][col][0];
        b = G.params->colorb[h][col][1];
        A = LT1B.params->colidx[a]; Asym = LT1B.params->qsym[a];
        B = R1A.params->colidx[b]; Bsym = R1A.params->qsym[b];
        if( ((Csym^Asym)==G_irr) && (Isym==Bsym))
          G.matrix[h][row][col] +=
            LT1B.matrix[Csym][C][A] * R1A.matrix[Isym][I][B];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);

  dpd_file2_mat_close(&LT1A);
  dpd_file2_mat_close(&LT1B);
  dpd_file2_close(&LT1A);
  dpd_file2_close(&LT1B);

  dpd_file2_mat_close(&R1A);
  dpd_file2_mat_close(&R1B);
  dpd_file2_close(&R1A);
  dpd_file2_close(&R1B);

  return;
}
