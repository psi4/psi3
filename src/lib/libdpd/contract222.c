#include <stdio.h>
#include <math.h>
#include <qt.h>
#include "dpd.h"

int dpd_contract222(dpdfile2 *X, dpdfile2 *Y, dpdfile2 *Z, int target_X,
		    int target_Y, double alpha, double beta)
{
  int h, nirreps, Xtrans, Ytrans, *numlinks;
#ifdef DPD_DEBUG
  int *xrow, *xcol, *yrow, *ycol, *zrow, *zcol;
#endif

  nirreps = X->params->nirreps;

  dpd_file2_mat_init(X);
  dpd_file2_mat_rd(X);
  dpd_file2_mat_init(Y);
  dpd_file2_mat_rd(Y);
  dpd_file2_mat_init(Z);
  if(fabs(beta) > 0.0) dpd_file2_mat_rd(Z);

  if(target_X == 0) { Xtrans = 0; numlinks = X->params->coltot; }
  else if(target_X == 1) { Xtrans = 1; numlinks = X->params->rowtot; }
  else {
      fprintf(stderr, "Junk X index %d in contract222\n", target_X);
      exit(target_X);
    }
  if(target_Y == 0) Ytrans = 1;
  else if(target_Y == 1) Ytrans = 0;
  else {
      fprintf(stderr, "Junk Y index %d in contract222\n", target_Y);
      exit(target_Y);
    }

#ifdef DPD_DEBUG
  if(Xtrans) { xrow = X->params->coltot; xcol = X->params->rowtot; }
  else { xrow = X->params->rowtot; xcol = X->params->coltot; }

  if(Ytrans) { yrow = Y->params->coltot; ycol = Y->params->rowtot; }
  else { yrow = Y->params->rowtot; ycol = Y->params->coltot; }

  zrow = Z->params->rowtot; zcol = Z->params->coltot;
  
  if((zrow != xrow) || (zcol != ycol) || (xcol != yrow)) {
      fprintf(stderr, "** Alignment error in contract222 **\n");
      dpd_error("dpd_contract222", stderr);
    }
#endif  

  for(h=0; h < nirreps; h++) {

      newmm(X->matrix[h], Xtrans, Y->matrix[h], Ytrans,
	    Z->matrix[h], Z->params->rowtot[h],
	    numlinks[h], Z->params->coltot[h], alpha, beta);
    }

  dpd_file2_mat_wrt(Z);
  dpd_file2_mat_close(X);
  dpd_file2_mat_close(Y);
  dpd_file2_mat_close(Z);

  return 0;
}
