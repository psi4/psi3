#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void taut_build(void)
{
  int h, ij, ab, i, j, a, b, I, J, A, B;
  int Isym, Jsym, Asym, Bsym;
  int nirreps;
  dpdbuf4 tautIJAB, tautijab, tautIjAb;
  dpdbuf4 tIJAB, tijab, tIjAb;
  dpdfile2 tIA, tia;

  nirreps = moinfo.nirreps;

  dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
  dpd_buf4_copy(&tIJAB, CC_TAMPS, "tautIJAB");
  dpd_buf4_close(&tIJAB);

  dpd_buf4_init(&tijab, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tijab");
  dpd_buf4_copy(&tijab, CC_TAMPS, "tautijab");
  dpd_buf4_close(&tijab);

  dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_copy(&tIjAb, CC_TAMPS, "tautIjAb");
  dpd_buf4_close(&tIjAb);

  dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_mat_init(&tIA);
  dpd_file2_mat_rd(&tIA);
  dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");
  dpd_file2_mat_init(&tia);
  dpd_file2_mat_rd(&tia);

  dpd_buf4_init(&tautIJAB, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tautIJAB");

  for(h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(&tautIJAB, h);
      dpd_buf4_mat_irrep_rd(&tautIJAB, h);

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

      dpd_buf4_mat_irrep_wrt(&tautIJAB, h);
      dpd_buf4_mat_irrep_close(&tautIJAB, h);
    }

  dpd_buf4_close(&tautIJAB);

  dpd_buf4_init(&tautijab, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tautijab");

  for(h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(&tautijab, h);
      dpd_buf4_mat_irrep_rd(&tautijab, h);

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

      dpd_buf4_mat_irrep_wrt(&tautijab, h);
      dpd_buf4_mat_irrep_close(&tautijab, h);
    }

  dpd_buf4_close(&tautijab);

  dpd_buf4_init(&tautIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tautIjAb");

  for(h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(&tautIjAb, h);
      dpd_buf4_mat_irrep_rd(&tautIjAb, h);

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

      dpd_buf4_mat_irrep_wrt(&tautIjAb, h);
      dpd_buf4_mat_irrep_close(&tautIjAb, h);
    }

  dpd_buf4_close(&tautIjAb);

  dpd_file2_mat_close(&tIA);
  dpd_file2_close(&tIA);
  dpd_file2_mat_close(&tia);
  dpd_file2_close(&tia);
}
