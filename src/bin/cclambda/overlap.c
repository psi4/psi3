#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void overlap(void)
{
  int h, nirreps;
  int row, col;
  int i,j,a,b,I,J,A,B,Isym,Jsym,Asym,Bsym;
  struct oe_dpdfile T1, L1, T1A, T1B;
  struct dpdbuf T2, L2;
  double value = 1.0;
  double ST1A, ST1B, ST2AA, ST2BB, ST2AB, ST12AA, ST12BB, ST12AB;

  nirreps = moinfo.nirreps;

  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "LIA", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  ST1A = dpd_oe_dot(&T1, &L1, 0, outfile);
  dpd_oe_file_close(&L1);
  dpd_oe_file_close(&T1);

  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "Lia", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  ST1B = dpd_oe_dot(&T1, &L1, 0, outfile);
  dpd_oe_file_close(&L1);
  dpd_oe_file_close(&T1);

  dpd_buf_init(&L2, CC_LAMPS, 2, 7, 2, 7, 0, "LIJAB", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 2, 7, 2, 7, 0, "tIJAB", 0, outfile);
  ST2AA = dpd_dot(&L2, &T2, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&L2);

  dpd_buf_init(&L2, CC_LAMPS, 2, 7, 2, 7, 0, "Lijab", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 2, 7, 2, 7, 0, "tijab", 0, outfile);
  ST2BB = dpd_dot(&L2, &T2, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&L2);

  dpd_buf_init(&L2, CC_LAMPS, 0, 5, 0, 5, 0, "LIjAb", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 0, 5, 0, 5, 0, "tIjAb", 0, outfile);
  ST2AB = dpd_dot(&L2, &T2, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&L2);

  dpd_oe_file_init(&T1A, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_mat_init(&T1A);
  dpd_oe_file_mat_rd(&T1A, 0, outfile);
  dpd_oe_file_init(&T1B, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_oe_file_mat_init(&T1B);
  dpd_oe_file_mat_rd(&T1B, 0, outfile);

  ST12AA = 0.0;
  dpd_buf_init(&L2, CC_LAMPS, 2, 7, 2, 7, 0, "LIJAB", 0, outfile);
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&L2, h);
      dpd_buf_mat_irrep_rd(&L2, h, 0, outfile);
      for(row=0; row < L2.params->rowtot[h]; row++) {
          i = L2.params->roworb[h][row][0];
          j = L2.params->roworb[h][row][1];
          I = T1A.params->rowidx[i]; Isym = T1A.params->psym[i];
          J = T1A.params->rowidx[j]; Jsym = T1A.params->psym[j];
          for(col=0; col < L2.params->coltot[h]; col++) {
              a = L2.params->colorb[h][col][0];
              b = L2.params->colorb[h][col][1];
              A = T1A.params->colidx[a]; Asym = T1A.params->qsym[a];
              B = T1A.params->colidx[b]; Bsym = T1A.params->qsym[b];
              if((Isym == Asym) && (Jsym == Bsym))
                 ST12AA += L2.matrix[h][row][col] * 
                          T1A.matrix[Isym][I][A] * T1A.matrix[Jsym][J][B];
              if((Isym == Bsym) && (Jsym == Asym))
                 ST12AA -= L2.matrix[h][row][col] * 
                          T1A.matrix[Isym][I][B] * T1A.matrix[Jsym][J][A];
            }
        }
      dpd_buf_mat_irrep_close(&L2, h);
    }
  dpd_buf_close(&L2);

  ST12BB = 0.0;

  dpd_buf_init(&L2, CC_LAMPS, 2, 7, 2, 7, 0, "Lijab", 0, outfile);
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&L2, h);
      dpd_buf_mat_irrep_rd(&L2, h, 0, outfile);
      for(row=0; row < L2.params->rowtot[h]; row++) {
          i = L2.params->roworb[h][row][0];
          j = L2.params->roworb[h][row][1];
          I = T1B.params->rowidx[i]; Isym = T1B.params->psym[i];
          J = T1B.params->rowidx[j]; Jsym = T1B.params->psym[j];
          for(col=0; col < L2.params->coltot[h]; col++) {
              a = L2.params->colorb[h][col][0];
              b = L2.params->colorb[h][col][1];
              A = T1B.params->colidx[a]; Asym = T1B.params->qsym[a];
              B = T1B.params->colidx[b]; Bsym = T1B.params->qsym[b];
              if((Isym == Asym) && (Jsym == Bsym))
                 ST12BB += L2.matrix[h][row][col] *
                          T1B.matrix[Isym][I][A] * T1B.matrix[Jsym][J][B];
              if((Isym == Bsym) && (Jsym == Asym))
                 ST12BB -= L2.matrix[h][row][col] *
                          T1B.matrix[Isym][I][B] * T1B.matrix[Jsym][J][A];
            }
        }
      dpd_buf_mat_irrep_close(&L2, h);
    }
  dpd_buf_close(&L2);

  ST12AB = 0.0;
  dpd_buf_init(&L2, CC_LAMPS, 0, 5, 0, 5, 0, "LIjAb", 0, outfile);
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&L2, h);
      dpd_buf_mat_irrep_rd(&L2, h, 0, outfile);
      for(row=0; row < L2.params->rowtot[h]; row++) {
          i = L2.params->roworb[h][row][0];
          j = L2.params->roworb[h][row][1];
          I = T1A.params->rowidx[i]; Isym = T1A.params->psym[i];
          J = T1B.params->rowidx[j]; Jsym = T1B.params->psym[j];
          for(col=0; col < L2.params->coltot[h]; col++) {
              a = L2.params->colorb[h][col][0];
              b = L2.params->colorb[h][col][1];
              A = T1A.params->colidx[a]; Asym = T1A.params->qsym[a];
              B = T1B.params->colidx[b]; Bsym = T1B.params->qsym[b];
              if((Isym == Asym) && (Jsym == Bsym))
                 ST12AB += L2.matrix[h][row][col] *
                          T1A.matrix[Isym][I][A] * T1B.matrix[Jsym][J][B];
            }
        }
      dpd_buf_mat_irrep_close(&L2, h);
    }
  dpd_buf_close(&L2);


  dpd_oe_file_mat_close(&T1A);
  dpd_oe_file_close(&T1A);
  dpd_oe_file_mat_close(&T1B);
  dpd_oe_file_close(&T1B);

/*
  fprintf(outfile, "\tST1A = %20.15f\n", ST1A);
  fprintf(outfile, "\tST1B = %20.15f\n", ST1B);
  fprintf(outfile, "\tST2AA = %20.15f\n", ST2AA);
  fprintf(outfile, "\tST2BB = %20.15f\n", ST2BB);
  fprintf(outfile, "\tST2AB = %20.15f\n", ST2AB);
  fprintf(outfile, "\tST12AA = %20.15f\n", ST12AA);
  fprintf(outfile, "\tST12BB = %20.15f\n", ST12BB);
  fprintf(outfile, "\tST12AB = %20.15f\n", ST12AB);
*/

  value = 1.0 - ST1A - ST1B - ST2AA - ST2BB - ST2AB + ST12AA + ST12BB + ST12AB;

  fprintf(outfile, "\tOverlap = %20.15f\n", value);
}
