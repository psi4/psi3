#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void tau_build(void)
{
  int h, ij, ab, i, j, a, b, I, J, A, B;
  int Isym, Jsym, Asym, Bsym;
  int nirreps;
  struct dpdbuf tauIJAB, tauijab, tauIjAb, tauiJaB, tauIjbA;
  struct dpdbuf tIJAB, tijab, tIjAb;
  struct oe_dpdfile tIA, tia;

  nirreps = moinfo.nirreps;

  dpd_buf_init(&tIJAB, CC_TAMPS, 2, 7, 2, 7, 0, "tIJAB", 0, outfile);
  dpd_copy(&tIJAB, CC_TAMPS, "tauIJAB", 0, outfile);
  dpd_buf_close(&tIJAB);

  dpd_buf_init(&tijab, CC_TAMPS, 2, 7, 2, 7, 0, "tijab", 0, outfile);
  dpd_copy(&tijab, CC_TAMPS, "tauijab", 0, outfile);
  dpd_buf_close(&tijab);

  dpd_buf_init(&tIjAb, CC_TAMPS, 0, 5, 0, 5, 0, "tIjAb", 0, outfile);
  dpd_copy(&tIjAb, CC_TAMPS, "tauIjAb", 0, outfile);
  dpd_buf_close(&tIjAb);

  dpd_oe_file_init(&tIA, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_mat_init(&tIA);
  dpd_oe_file_mat_rd(&tIA, 0, outfile);
  dpd_oe_file_init(&tia, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_oe_file_mat_init(&tia);
  dpd_oe_file_mat_rd(&tia, 0, outfile);

  dpd_buf_init(&tauIJAB, CC_TAMPS, 2, 7, 2, 7, 0, "tauIJAB", 0, outfile);

  for(h=0; h < nirreps; h++) {

      dpd_buf_mat_irrep_init(&tauIJAB, h);
      dpd_buf_mat_irrep_rd(&tauIJAB, h, 0, outfile);

      for(ij=0; ij < tauIJAB.params->rowtot[h]; ij++) {
	  i = tauIJAB.params->roworb[h][ij][0];
	  j = tauIJAB.params->roworb[h][ij][1];
	  I = tIA.params->rowidx[i];
	  J = tIA.params->rowidx[j];
	  Isym = tIA.params->psym[i];
	  Jsym = tIA.params->psym[j];
	  for(ab=0; ab < tauIJAB.params->coltot[h]; ab++) {
	      a = tauIJAB.params->colorb[h][ab][0];
	      b = tauIJAB.params->colorb[h][ab][1];
	      A = tIA.params->colidx[a];
	      B = tIA.params->colidx[b];
	      Asym = tIA.params->qsym[a];
	      Bsym = tIA.params->qsym[b];

	      if((Isym==Asym) && (Jsym==Bsym))
		  tauIJAB.matrix[h][ij][ab] +=
		      (tIA.matrix[Isym][I][A] * tIA.matrix[Jsym][J][B]);
	      if((Isym==Bsym) && (Jsym==Asym))
		  tauIJAB.matrix[h][ij][ab] -=
		      (tIA.matrix[Isym][I][B] * tIA.matrix[Jsym][J][A]);

	    }
	}

      dpd_buf_mat_irrep_wrt(&tauIJAB, h, 0, outfile);
      dpd_buf_mat_irrep_close(&tauIJAB, h);
    }

  dpd_buf_close(&tauIJAB);

  dpd_buf_init(&tauijab, CC_TAMPS, 2, 7, 2, 7, 0, "tauijab", 0, outfile);

  for(h=0; h < nirreps; h++) {

      dpd_buf_mat_irrep_init(&tauijab, h);
      dpd_buf_mat_irrep_rd(&tauijab, h, 0, outfile);

      for(ij=0; ij < tauijab.params->rowtot[h]; ij++) {
	  i = tauijab.params->roworb[h][ij][0];
	  j = tauijab.params->roworb[h][ij][1];
	  I = tia.params->rowidx[i];
	  J = tia.params->rowidx[j];
	  Isym = tia.params->psym[i];
	  Jsym = tia.params->psym[j];
	  for(ab=0; ab < tauijab.params->coltot[h]; ab++) {
	      a = tauijab.params->colorb[h][ab][0];
	      b = tauijab.params->colorb[h][ab][1];
	      A = tia.params->colidx[a];
	      B = tia.params->colidx[b];
	      Asym = tia.params->qsym[a];
	      Bsym = tia.params->qsym[b];

	      if((Isym==Asym) && (Jsym==Bsym))
		  tauijab.matrix[h][ij][ab] +=
		      (tia.matrix[Isym][I][A] * tia.matrix[Jsym][J][B]);
	      if((Isym==Bsym) && (Jsym==Asym))
		  tauijab.matrix[h][ij][ab] -=
		      (tia.matrix[Isym][I][B] * tia.matrix[Jsym][J][A]);

	    }
	}

      dpd_buf_mat_irrep_wrt(&tauijab, h, 0, outfile);
      dpd_buf_mat_irrep_close(&tauijab, h);
    }

  dpd_buf_close(&tauijab);

  dpd_buf_init(&tauIjAb, CC_TAMPS, 0, 5, 0, 5, 0, "tauIjAb", 0, outfile);

  for(h=0; h < nirreps; h++) {

      dpd_buf_mat_irrep_init(&tauIjAb, h);
      dpd_buf_mat_irrep_rd(&tauIjAb, h, 0, outfile);

      for(ij=0; ij < tauIjAb.params->rowtot[h]; ij++) {
	  i = tauIjAb.params->roworb[h][ij][0];
	  j = tauIjAb.params->roworb[h][ij][1];
	  I = tIA.params->rowidx[i];
	  J = tia.params->rowidx[j];
	  Isym = tIA.params->psym[i];
	  Jsym = tia.params->psym[j];
	  for(ab=0; ab < tauIjAb.params->coltot[h]; ab++) {
	      a = tauIjAb.params->colorb[h][ab][0];
	      b = tauIjAb.params->colorb[h][ab][1];
	      A = tIA.params->colidx[a];
	      B = tia.params->colidx[b];
	      Asym = tIA.params->qsym[a];
	      Bsym = tia.params->qsym[b];

	      if((Isym==Asym) && (Jsym==Bsym))
		  tauIjAb.matrix[h][ij][ab] +=
		      (tIA.matrix[Isym][I][A] * tia.matrix[Jsym][J][B]);

	    }
	}

      dpd_buf_mat_irrep_wrt(&tauIjAb, h, 0, outfile);
      dpd_buf_mat_irrep_close(&tauIjAb, h);
    }

  /* This will generate the tauBA and tauIjbA files from tauIjAb */
  dpd_swap34(&tauIjAb, CC_TAMPS, 0, 5, "tauIjbA", 0, outfile);
  dpd_buf_close(&tauIjAb);
  dpd_buf_init(&tauIjbA, CC_TAMPS, 0, 5, 0, 5, 0, "tauIjbA", 0, outfile);
  dpd_swap12(&tauIjbA, CC_TAMPS, 0, 5, "tauiJaB", 0, outfile);
  dpd_buf_close(&tauIjbA);

  dpd_oe_file_mat_close(&tIA);
  dpd_oe_file_close(&tIA);
  dpd_oe_file_mat_close(&tia);
  dpd_oe_file_close(&tia);
}
