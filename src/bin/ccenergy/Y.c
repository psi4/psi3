#include <stdio.h>
#include <dpd.h>
#include <qt.h>
#define EXTERN
#include "globals.h"

void Y_build(void)
{
  int h, nirreps;
  int ia, jb, i, a, j, b, I, A, J, B, Isym, Asym, Jsym, Bsym;
  struct dpdbuf Y;
  struct oe_dpdfile tIA, tia;

/*  timer_on("Y_build", outfile); */

  nirreps = moinfo.nirreps;

  dpd_oe_file_init(&tIA, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_mat_init(&tIA);
  dpd_oe_file_mat_rd(&tIA, 0, outfile);
  dpd_oe_file_init(&tia, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_oe_file_mat_init(&tia);
  dpd_oe_file_mat_rd(&tia, 0, outfile);

  dpd_buf_init(&Y, CC_MISC, 10, 10, 10, 10, 0, "YIAJB", 0, outfile);

  for(h=0; h < nirreps; h++) {

      dpd_buf_mat_irrep_init(&Y, h);

      for(ia=0; ia < Y.params->rowtot[h]; ia++) {
	  i = Y.params->roworb[h][ia][0];
	  a = Y.params->roworb[h][ia][1];
	  I = tIA.params->rowidx[i];  Isym = tIA.params->psym[i];
	  A = tIA.params->colidx[a];  Asym = tIA.params->qsym[a];
	  for(jb=0; jb < Y.params->coltot[h]; jb++) {
	      j = Y.params->colorb[h][jb][0];
	      b = Y.params->colorb[h][jb][1];
	      J = tIA.params->rowidx[j];  Jsym = tIA.params->psym[j];
	      B = tIA.params->colidx[b];  Bsym = tIA.params->qsym[b];

	      if((Isym == Bsym) && (Jsym == Asym))
		  Y.matrix[h][ia][jb] =
			tIA.matrix[Isym][I][B] * tIA.matrix[Jsym][J][A];
	      
	    }
	}

      dpd_buf_mat_irrep_wrt(&Y, h, 0, outfile);
      dpd_buf_mat_irrep_close(&Y, h);
      
    }

  dpd_buf_close(&Y);

  dpd_buf_init(&Y, CC_MISC, 10, 10, 10, 10, 0, "Yiajb", 0, outfile);

  for(h=0; h < nirreps; h++) {

      dpd_buf_mat_irrep_init(&Y, h);

      for(ia=0; ia < Y.params->rowtot[h]; ia++) {
	  i = Y.params->roworb[h][ia][0];
	  a = Y.params->roworb[h][ia][1];
	  I = tia.params->rowidx[i];  Isym = tia.params->psym[i];
	  A = tia.params->colidx[a];  Asym = tia.params->qsym[a];
	  for(jb=0; jb < Y.params->coltot[h]; jb++) {
	      j = Y.params->colorb[h][jb][0];
	      b = Y.params->colorb[h][jb][1];
	      J = tia.params->rowidx[j];  Jsym = tia.params->psym[j];
	      B = tia.params->colidx[b];  Bsym = tia.params->qsym[b];

	      if((Isym == Bsym) && (Jsym == Asym))
		  Y.matrix[h][ia][jb] =
		      tia.matrix[Isym][I][B] * tia.matrix[Jsym][J][A];
	      
	    }
	}

      dpd_buf_mat_irrep_wrt(&Y, h, 0, outfile);
      dpd_buf_mat_irrep_close(&Y, h);
      
    }

  dpd_buf_close(&Y);

  dpd_buf_init(&Y, CC_MISC, 10, 10, 10, 10, 0, "YIajB", 0, outfile);

  for(h=0; h < nirreps; h++) {

      dpd_buf_mat_irrep_init(&Y, h);

      for(ia=0; ia < Y.params->rowtot[h]; ia++) {
	  i = Y.params->roworb[h][ia][0];
	  a = Y.params->roworb[h][ia][1];
	  I = tIA.params->rowidx[i];  Isym = tIA.params->psym[i];
	  A = tia.params->colidx[a];  Asym = tia.params->qsym[a];
	  for(jb=0; jb < Y.params->coltot[h]; jb++) {
	      j = Y.params->colorb[h][jb][0];
	      b = Y.params->colorb[h][jb][1];
	      J = tia.params->rowidx[j];  Jsym = tia.params->psym[j];
	      B = tIA.params->colidx[b];  Bsym = tIA.params->qsym[b];

	      if((Isym == Bsym) && (Jsym == Asym))
		  Y.matrix[h][ia][jb] =
			tIA.matrix[Isym][I][B] * tia.matrix[Jsym][J][A];
	      
	    }
	}

      dpd_buf_mat_irrep_wrt(&Y, h, 0, outfile);
      dpd_buf_mat_irrep_close(&Y, h);
      
    }

  dpd_buf_close(&Y);

  dpd_buf_init(&Y, CC_MISC, 10, 10, 10, 10, 0, "YiAJb", 0, outfile);

  for(h=0; h < nirreps; h++) {

      dpd_buf_mat_irrep_init(&Y, h);

      for(ia=0; ia < Y.params->rowtot[h]; ia++) {
	  i = Y.params->roworb[h][ia][0];
	  a = Y.params->roworb[h][ia][1];
	  I = tia.params->rowidx[i];  Isym = tia.params->psym[i];
	  A = tIA.params->colidx[a];  Asym = tIA.params->qsym[a];
	  for(jb=0; jb < Y.params->coltot[h]; jb++) {
	      j = Y.params->colorb[h][jb][0];
	      b = Y.params->colorb[h][jb][1];
	      J = tIA.params->rowidx[j];  Jsym = tIA.params->psym[j];
	      B = tia.params->colidx[b];  Bsym = tia.params->qsym[b];

	      if((Isym == Bsym) && (Jsym == Asym))
		  Y.matrix[h][ia][jb] =
			tia.matrix[Isym][I][B] * tIA.matrix[Jsym][J][A];
	      
	    }
	}

      dpd_buf_mat_irrep_wrt(&Y, h, 0, outfile);
      dpd_buf_mat_irrep_close(&Y, h);
      
    }

  dpd_buf_close(&Y);

  dpd_oe_file_mat_close(&tIA);
  dpd_oe_file_close(&tIA);
  dpd_oe_file_mat_close(&tia);
  dpd_oe_file_close(&tia);

/*  timer_off("Y_build", outfile); */
}
