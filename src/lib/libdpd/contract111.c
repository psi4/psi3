#include <stdio.h>
#include <math.h>
#include <qt.h>
#include "dpd.h"

int dpd_contract111(struct oe_dpdfile *X, struct oe_dpdfile *Y, struct
		    oe_dpdfile *Z, int target_X, int target_Y, double
		    alpha, double beta, int print_flag, FILE *outfile)
{
  int h, nirreps, Xtrans, Ytrans, *numlinks;
#ifdef DPD_DEBUG
  int *xrow, *xcol, *yrow, *ycol, *zrow, *zcol;
#endif

  nirreps = X->params->nirreps;

  dpd_oe_file_mat_init(X);
  dpd_oe_file_mat_rd(X, print_flag, outfile);
  dpd_oe_file_mat_init(Y);
  dpd_oe_file_mat_rd(Y, print_flag, outfile);
  dpd_oe_file_mat_init(Z);
  if(fabs(beta) > 0.0) dpd_oe_file_mat_rd(Z, print_flag, outfile);

  if(target_X == 0) { Xtrans = 0; numlinks = X->params->coltot; }
  else if(target_X == 1) { Xtrans = 1; numlinks = X->params->rowtot; }
  else {
      fprintf(outfile, "Junk X index %d in contract111\n", target_X);
      exit(target_X);
    }
  if(target_Y == 0) Ytrans = 1;
  else if(target_Y == 1) Ytrans = 0;
  else {
      fprintf(outfile, "Junk Y index %d in contract111\n", target_Y);
      exit(target_Y);
    }

#ifdef DPD_DEBUG
  if(Xtrans) { xrow = X->params->coltot; xcol = X->params->rowtot; }
  else { xrow = X->params->rowtot; xcol = X->params->coltot; }

  if(Ytrans) { yrow = Y->params->coltot; ycol = Y->params->rowtot; }
  else { yrow = Y->params->rowtot; ycol = Y->params->coltot; }

  zrow = Z->params->rowtot; zcol = Z->params->coltot;
  
  if((zrow != xrow) || (zcol != ycol) || (xcol != yrow)) {
      fprintf(outfile, "** Alignment error in contract111 **\n");
      dpd_error("dpd_contract111", outfile);
    }
#endif  

  for(h=0; h < nirreps; h++) {

      newmm(X->matrix[h], Xtrans, Y->matrix[h], Ytrans,
	    Z->matrix[h], Z->params->rowtot[h],
	    numlinks[h], Z->params->coltot[h], alpha, beta);
    }

  dpd_oe_file_mat_wrt(Z, print_flag, outfile);
  dpd_oe_file_mat_close(X);
  dpd_oe_file_mat_close(Y);
  dpd_oe_file_mat_close(Z);

  return 0;
}
