#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void taut_build(void)
{
  int h, ij, ab, i, j, a, b, I, J, A, B;
  int Isym, Jsym, Asym, Bsym;
  int nirreps;
  struct dpdbuf tautIJAB, tautijab, tautIjAb;
  struct dpdbuf tIJAB, tijab, tIjAb;
  struct oe_dpdfile tIA, tia;

  nirreps = moinfo.nirreps;

  dpd_buf_init(&tIJAB, CC_TAMPS, 2, 7, 2, 7, 0, "tIJAB", 0, outfile);
  dpd_copy(&tIJAB, CC_TAMPS, "tautIJAB", 0, outfile);
  dpd_buf_close(&tIJAB);

  dpd_buf_init(&tijab, CC_TAMPS, 2, 7, 2, 7, 0, "tijab", 0, outfile);
  dpd_copy(&tijab, CC_TAMPS, "tautijab", 0, outfile);
  dpd_buf_close(&tijab);

  dpd_buf_init(&tIjAb, CC_TAMPS, 0, 5, 0, 5, 0, "tIjAb", 0, outfile);
  dpd_copy(&tIjAb, CC_TAMPS, "tautIjAb", 0, outfile);
  dpd_buf_close(&tIjAb);

  dpd_oe_file_init(&tIA, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_mat_init(&tIA);
  dpd_oe_file_mat_rd(&tIA, 0, outfile);
  dpd_oe_file_init(&tia, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_oe_file_mat_init(&tia);
  dpd_oe_file_mat_rd(&tia, 0, outfile);

  dpd_buf_init(&tautIJAB, CC_TAMPS, 2, 7, 2, 7, 0, "tautIJAB", 0, outfile);

  for(h=0; h < nirreps; h++) {

      dpd_buf_mat_irrep_init(&tautIJAB, h);
      dpd_buf_mat_irrep_rd(&tautIJAB, h, 0, outfile);

      for(ij=0; ij < tautIJAB.params->rowtot[h]; ij++) {
	  i = tautIJAB.params->roworb[h][ij][0];
	  j = tautIJAB.params->roworb[h][ij][1];
	  I = tIA.params->rowidx[i];
	  J = tIA.params->rowidx[j];
	  Isym = tIA.params->psym[i];
	  Jsym = tIA.params->psym[j];
	  for(ab=0; ab < tautIJAB.params->coltot[h]; ab++) {
	      a = tautIJAB.params->colorb[h][ab][0];
	      b = tautIJAB.params->colorb[h][ab][1];
	      A = tIA.params->colidx[a];
	      B = tIA.params->colidx[b];
	      Asym = tIA.params->qsym[a];
	      Bsym = tIA.params->qsym[b];

	      if((Isym==Asym) && (Jsym==Bsym))
		  tautIJAB.matrix[h][ij][ab] +=
		   0.5*(tIA.matrix[Isym][I][A] * tIA.matrix[Jsym][J][B]);
	      if((Isym==Bsym) && (Jsym==Asym))
		  tautIJAB.matrix[h][ij][ab] -=
		   0.5*(tIA.matrix[Isym][I][B] * tIA.matrix[Jsym][J][A]);

	    }
	}

      dpd_buf_mat_irrep_wrt(&tautIJAB, h, 0, outfile);
      dpd_buf_mat_irrep_close(&tautIJAB, h);
    }

  dpd_buf_close(&tautIJAB);

  dpd_buf_init(&tautijab, CC_TAMPS, 2, 7, 2, 7, 0, "tautijab", 0, outfile);

  for(h=0; h < nirreps; h++) {

      dpd_buf_mat_irrep_init(&tautijab, h);
      dpd_buf_mat_irrep_rd(&tautijab, h, 0, outfile);

      for(ij=0; ij < tautijab.params->rowtot[h]; ij++) {
	  i = tautijab.params->roworb[h][ij][0];
	  j = tautijab.params->roworb[h][ij][1];
	  I = tia.params->rowidx[i];
	  J = tia.params->rowidx[j];
	  Isym = tia.params->psym[i];
	  Jsym = tia.params->psym[j];
	  for(ab=0; ab < tautijab.params->coltot[h]; ab++) {
	      a = tautijab.params->colorb[h][ab][0];
	      b = tautijab.params->colorb[h][ab][1];
	      A = tia.params->colidx[a];
	      B = tia.params->colidx[b];
	      Asym = tia.params->qsym[a];
	      Bsym = tia.params->qsym[b];

	      if((Isym==Asym) && (Jsym==Bsym))
		  tautijab.matrix[h][ij][ab] +=
		   0.5*(tia.matrix[Isym][I][A] * tia.matrix[Jsym][J][B]);
	      if((Isym==Bsym) && (Jsym==Asym))
		  tautijab.matrix[h][ij][ab] -=
		   0.5*(tia.matrix[Isym][I][B] * tia.matrix[Jsym][J][A]);

	    }
	}

      dpd_buf_mat_irrep_wrt(&tautijab, h, 0, outfile);
      dpd_buf_mat_irrep_close(&tautijab, h);
    }

  dpd_buf_close(&tautijab);

  dpd_buf_init(&tautIjAb, CC_TAMPS, 0, 5, 0, 5, 0, "tautIjAb", 0, outfile);

  for(h=0; h < nirreps; h++) {

      dpd_buf_mat_irrep_init(&tautIjAb, h);
      dpd_buf_mat_irrep_rd(&tautIjAb, h, 0, outfile);

      for(ij=0; ij < tautIjAb.params->rowtot[h]; ij++) {
	  i = tautIjAb.params->roworb[h][ij][0];
	  j = tautIjAb.params->roworb[h][ij][1];
	  I = tIA.params->rowidx[i];
	  J = tia.params->rowidx[j];
	  Isym = tIA.params->psym[i];
	  Jsym = tia.params->psym[j];
	  for(ab=0; ab < tautIjAb.params->coltot[h]; ab++) {
	      a = tautIjAb.params->colorb[h][ab][0];
	      b = tautIjAb.params->colorb[h][ab][1];
	      A = tIA.params->colidx[a];
	      B = tia.params->colidx[b];
	      Asym = tIA.params->qsym[a];
	      Bsym = tia.params->qsym[b];

	      if((Isym==Asym) && (Jsym==Bsym))
		  tautIjAb.matrix[h][ij][ab] +=
		   0.5*(tIA.matrix[Isym][I][A] * tia.matrix[Jsym][J][B]);

	    }
	}

      dpd_buf_mat_irrep_wrt(&tautIjAb, h, 0, outfile);
      dpd_buf_mat_irrep_close(&tautIjAb, h);
    }

  dpd_buf_close(&tautIjAb);

  dpd_oe_file_mat_close(&tIA);
  dpd_oe_file_close(&tIA);
  dpd_oe_file_mat_close(&tia);
  dpd_oe_file_close(&tia);
}
