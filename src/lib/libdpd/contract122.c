#include <stdio.h>
#include <math.h>
#include <qt.h>
#include "dpd.h"

/*
** target_X  = which index on X is external 
** target_Y  = which index on Y is external
**
** Note that the external index on X MUST be the left index in Z in this version.
*/

int dpd_contract122(struct dpdbuf *X, struct dpdbuf *Y, struct oe_dpdfile *Z,
		    int target_X, int target_Y, double alpha, double beta,
		    int print_flag, FILE *outfile)
{
  int h,nirreps,Gtar;
  int rking=0;
  struct dpdtrans Xt, Yt;
  double ***Xmat, ***Ymat, ***Zmat;
  int Xtrans, Ytrans, *numlinks;
#ifdef DPD_DEBUG
  int *xrow, *xcol, *yrow, *ycol, *zrow, *zcol;
#endif

  nirreps = X->params->nirreps;

  if((target_X == 1) || (target_X == 2)) dpd_trans_init(&Xt, X);

  if((target_Y == 1) || (target_Y == 2)) dpd_trans_init(&Yt, Y);

  if(fabs(beta) > 0.0) dpd_oe_scm(Z, beta, print_flag, outfile);
  dpd_oe_file_mat_init(Z);
  if(fabs(beta) > 0.0) dpd_oe_file_mat_rd(Z, print_flag, outfile);

#ifdef DPD_DEBUG
  zrow = Z->params->rowtot;
  zcol = Z->params->coltot;
#endif

  for(h=0; h < nirreps; h++) {
      if(target_X == 0) {
	  dpd_buf_mat_irrep_init(X, h);
	  dpd_buf_mat_irrep_rd(X, h, print_flag, outfile);
	  dpd_buf_mat_irrep_shift13(X, h, print_flag, outfile);
	  Xmat = X->shift.matrix[h];
	  Xtrans = 0;
	  numlinks = X->shift.coltot[h];
#ifdef DPD_DEBUG
	  xrow = X->shift.rowtot[h];
	  xcol = X->shift.coltot[h];
#endif
	}
      else if(target_X == 1) {
	  dpd_buf_mat_irrep_init(X, h);
	  dpd_buf_mat_irrep_rd(X, h, print_flag, outfile);
	  dpd_trans_mat_irrep_init(&Xt, h);
	  dpd_trans_mat_irrep_rd(&Xt, h);
	  dpd_buf_mat_irrep_close(X, h);
	  dpd_trans_mat_irrep_shift31(&Xt, h, print_flag, outfile);
	  rking = 1;
	  Xmat = Xt.shift.matrix[h];
	  Xtrans = 1;
	  numlinks = Xt.shift.rowtot[h];
#ifdef DPD_DEBUG	  
	  xrow = Xt.shift.coltot[h];
	  xcol = Xt.shift.rowtot[h];
#endif
	}
      else if(target_X == 2) {
	  dpd_buf_mat_irrep_init(X, h);
	  dpd_buf_mat_irrep_rd(X, h, print_flag, outfile);
	  dpd_trans_mat_irrep_init(&Xt, h);
	  dpd_trans_mat_irrep_rd(&Xt, h);
	  dpd_buf_mat_irrep_close(X, h);
	  dpd_trans_mat_irrep_shift13(&Xt, h, print_flag, outfile);
	  Xmat = Xt.shift.matrix[h];
	  Xtrans = 0;
	  numlinks = Xt.shift.coltot[h];
#ifdef DPD_DEBUG	  
	  xrow = Xt.shift.rowtot[h];
	  xcol = Xt.shift.coltot[h];
#endif	  
	}
      else if(target_X == 3) {
	  dpd_buf_mat_irrep_init(X, h);
	  dpd_buf_mat_irrep_rd(X, h, print_flag, outfile);
	  dpd_buf_mat_irrep_shift31(X, h, print_flag, outfile);
	  rking = 1;
	  Xmat = X->shift.matrix[h];
	  Xtrans = 1;
	  numlinks = X->shift.rowtot[h];
#ifdef DPD_DEBUG	  
	  xrow = X->shift.coltot[h];
	  xcol = X->shift.rowtot[h];
#endif	    
	}
      else {
	  fprintf(outfile, "Junk X index %d in contract122\n", target_X);
	  exit(target_X);
	}
      if(target_Y == 0) {
	  dpd_buf_mat_irrep_init(Y, h);
	  dpd_buf_mat_irrep_rd(Y, h, print_flag, outfile);
	  dpd_buf_mat_irrep_shift13(Y, h, print_flag, outfile);
	  Ymat = Y->shift.matrix[h];
	  Ytrans = 1;
#ifdef DPD_DEBUG	  
	  yrow = Y->shift.coltot[h];
	  ycol = Y->shift.rowtot[h];
#endif		  
	}
      else if(target_Y == 1) {
	  dpd_buf_mat_irrep_init(Y, h);
	  dpd_buf_mat_irrep_rd(Y, h, print_flag, outfile);
	  dpd_trans_mat_irrep_init(&Yt, h);
	  dpd_trans_mat_irrep_rd(&Yt, h);
	  dpd_buf_mat_irrep_close(Y, h);
	  dpd_trans_mat_irrep_shift31(&Yt, h, print_flag, outfile);
	  rking = 1;
	  Ymat = Yt.shift.matrix[h];
	  Ytrans = 0;
#ifdef DPD_DEBUG	  
	  yrow = Yt.shift.rowtot[h];
	  ycol = Yt.shift.coltot[h];
#endif	  
	}
      else if(target_Y == 2) {
	  dpd_buf_mat_irrep_init(Y, h);
	  dpd_buf_mat_irrep_rd(Y, h, print_flag, outfile);
	  dpd_trans_mat_irrep_init(&Yt, h);
	  dpd_trans_mat_irrep_rd(&Yt, h);
	  dpd_buf_mat_irrep_close(Y, h);
	  dpd_trans_mat_irrep_shift13(&Yt, h, print_flag, outfile);
	  Ymat = Yt.shift.matrix[h];
	  Ytrans = 1;
#ifdef DPD_DEBUG	  
	  yrow = Yt.shift.coltot[h];
	  ycol = Yt.shift.rowtot[h];
#endif	
	}
      else if(target_Y == 3) {
	  dpd_buf_mat_irrep_init(Y, h);
	  dpd_buf_mat_irrep_rd(Y, h, print_flag, outfile);
	  dpd_buf_mat_irrep_shift31(Y, h, print_flag, outfile);
	  rking = 1;
	  Ymat = Y->shift.matrix[h];
	  Ytrans = 0;
#ifdef DPD_DEBUG	  
	  yrow = Y->shift.rowtot[h];
	  ycol = Y->shift.coltot[h];
#endif	  
	}
      else {
	  fprintf(outfile, "Junk Y index %d in contract122\n", target_Y);
	  exit(target_Y);
	}

      if(rking)
	  for(Gtar=0; Gtar < nirreps; Gtar++) {
#ifdef DPD_DEBUG
	      if((xrow[Gtar] != zrow[Gtar]) || (ycol[Gtar] != zcol[Gtar]) ||
		 (xcol[Gtar] != yrow[Gtar])) {
		  fprintf(outfile, "** Alignment error in contract122 **\n");
		  fprintf(outfile, "** Irrep %d; Subirrep %d **\n", h,Gtar);
		  dpd_error("dpd_contract122", outfile);
		}
#endif
	      newmm_rking(Xmat[Gtar], Xtrans, Ymat[Gtar], Ytrans,
			  Z->matrix[Gtar], Z->params->rowtot[Gtar],
			  numlinks[Gtar], Z->params->coltot[Gtar],
			  alpha, 1.0);
	    }
      else
	  for(Gtar=0; Gtar < nirreps; Gtar++) {
#ifdef DPD_DEBUG
	      if((xrow[Gtar] != zrow[Gtar]) || (ycol[Gtar] != zcol[Gtar]) ||
		 (xcol[Gtar] != yrow[Gtar])) {
		  fprintf(outfile, "** Alignment error in contract122 **\n");
		  fprintf(outfile, "** Irrep %d; Subirrep %d **\n", h,Gtar);
		  dpd_error("dpd_contract122", outfile);
		}
#endif	      
	      newmm(Xmat[Gtar], Xtrans, Ymat[Gtar], Ytrans,
		    Z->matrix[Gtar], Z->params->rowtot[Gtar],
		    numlinks[Gtar], Z->params->coltot[Gtar],
		    alpha, 1.0);
	    }


      if(target_X == 0) dpd_buf_mat_irrep_close(X, h);
      else if(target_X == 1) dpd_trans_mat_irrep_close(&Xt, h);
      else if(target_X == 2) dpd_trans_mat_irrep_close(&Xt, h);
      else if(target_X == 3) dpd_buf_mat_irrep_close(X, h);

      if(target_Y == 0) dpd_buf_mat_irrep_close(Y, h);
      else if(target_Y == 1) dpd_trans_mat_irrep_close(&Yt, h);
      else if(target_Y == 2) dpd_trans_mat_irrep_close(&Yt, h);
      else if(target_Y == 3) dpd_buf_mat_irrep_close(Y, h);
    }

  if((target_X == 1) || (target_X == 2)) dpd_trans_close(&Xt);
  if((target_Y == 1) || (target_Y == 2)) dpd_trans_close(&Yt);


  dpd_oe_file_mat_wrt(Z, print_flag, outfile);
  dpd_oe_file_mat_close(Z);

  return 0;
}


