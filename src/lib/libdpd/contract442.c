#include <stdio.h>
#include <math.h>
#include <libqt/qt.h>
#include "dpd.h"

/* dpd_contract442(): Contracts a four-index quantity with another
** four-index quantity to compute the contribution to a
** two-index quantity, beta * Z = alpha * X * Y.
**
** Arguments:
**   dpdbuf4 *X: A pointer to the left four-index file.
**   dpdbuf4 *Y: A pointer to the four-index buffer.
**   dpdfile2 *Z: A pointer to the target two-index buffer.
**   int target_X: Indicates which index on X is to the target (takes a value of
**              0, 1, 2, or 3).
**   int target_Y: Indicates which index on Y is to the target (takes a value of
**              0, 1, 2, or 3).
**   double alpha: A prefactor for the product alpha * X * Y.
**   double beta: A prefactor for the target beta * Z.
*/

int dpd_contract442(dpdbuf4 *X, dpdbuf4 *Y, dpdfile2 *Z, int target_X,
		    int target_Y, double alpha, double beta)
{
  int h,nirreps,Gtar;
  int rking=0;
  dpdtrans4 Xt, Yt;
  double ***Xmat, ***Ymat, ***Zmat;
  int Xtrans, Ytrans, *numlinks;
#ifdef DPD_DEBUG
  int *xrow, *xcol, *yrow, *ycol, *zrow, *zcol;
#endif

  nirreps = X->params->nirreps;

  if((target_X == 1) || (target_X == 2)) dpd_trans4_init(&Xt, X);

  if((target_Y == 1) || (target_Y == 2)) dpd_trans4_init(&Yt, Y);

/*  if(fabs(beta) > 0.0) dpd_file2_scm(Z, beta); */
  dpd_file2_scm(Z, beta);
  dpd_file2_mat_init(Z);
/*  if(fabs(beta) > 0.0) dpd_file2_mat_rd(Z); */
  dpd_file2_mat_rd(Z);

#ifdef DPD_DEBUG
  zrow = Z->params->rowtot;
  zcol = Z->params->coltot;
#endif

  for(h=0; h < nirreps; h++) {
      if(target_X == 0) {
	  dpd_buf4_mat_irrep_init(X, h);
	  dpd_buf4_mat_irrep_rd(X, h);
	  dpd_buf4_mat_irrep_shift13(X, h);
	  Xmat = X->shift.matrix[h];
	  Xtrans = 0;
	  numlinks = X->shift.coltot[h];
#ifdef DPD_DEBUG
	  xrow = X->shift.rowtot[h];
	  xcol = X->shift.coltot[h];
#endif
	}
      else if(target_X == 1) {
	  dpd_buf4_mat_irrep_init(X, h);
	  dpd_buf4_mat_irrep_rd(X, h);
	  dpd_trans4_mat_irrep_init(&Xt, h);
	  dpd_trans4_mat_irrep_rd(&Xt, h);
	  dpd_buf4_mat_irrep_close(X, h);
	  dpd_trans4_mat_irrep_shift31(&Xt, h);
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
	  dpd_buf4_mat_irrep_init(X, h);
	  dpd_buf4_mat_irrep_rd(X, h);
	  dpd_trans4_mat_irrep_init(&Xt, h);
	  dpd_trans4_mat_irrep_rd(&Xt, h);
	  dpd_buf4_mat_irrep_close(X, h);
	  dpd_trans4_mat_irrep_shift13(&Xt, h);
	  Xmat = Xt.shift.matrix[h];
	  Xtrans = 0;
	  numlinks = Xt.shift.coltot[h];
#ifdef DPD_DEBUG	  
	  xrow = Xt.shift.rowtot[h];
	  xcol = Xt.shift.coltot[h];
#endif	  
	}
      else if(target_X == 3) {
	  dpd_buf4_mat_irrep_init(X, h);
	  dpd_buf4_mat_irrep_rd(X, h);
	  dpd_buf4_mat_irrep_shift31(X, h);
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
	  fprintf(stderr, "Junk X index %d in dpd_contract442\n", target_X);
	  exit(target_X);
	}
      if(target_Y == 0) {
	  dpd_buf4_mat_irrep_init(Y, h);
	  dpd_buf4_mat_irrep_rd(Y, h);
	  dpd_buf4_mat_irrep_shift13(Y, h);
	  Ymat = Y->shift.matrix[h];
	  Ytrans = 1;
#ifdef DPD_DEBUG	  
	  yrow = Y->shift.coltot[h];
	  ycol = Y->shift.rowtot[h];
#endif		  
	}
      else if(target_Y == 1) {
	  dpd_buf4_mat_irrep_init(Y, h);
	  dpd_buf4_mat_irrep_rd(Y, h);
	  dpd_trans4_mat_irrep_init(&Yt, h);
	  dpd_trans4_mat_irrep_rd(&Yt, h);
	  dpd_buf4_mat_irrep_close(Y, h);
	  dpd_trans4_mat_irrep_shift31(&Yt, h);
	  rking = 1;
	  Ymat = Yt.shift.matrix[h];
	  Ytrans = 0;
#ifdef DPD_DEBUG	  
	  yrow = Yt.shift.rowtot[h];
	  ycol = Yt.shift.coltot[h];
#endif	  
	}
      else if(target_Y == 2) {
	  dpd_buf4_mat_irrep_init(Y, h);
	  dpd_buf4_mat_irrep_rd(Y, h);
	  dpd_trans4_mat_irrep_init(&Yt, h);
	  dpd_trans4_mat_irrep_rd(&Yt, h);
	  dpd_buf4_mat_irrep_close(Y, h);
	  dpd_trans4_mat_irrep_shift13(&Yt, h);
	  Ymat = Yt.shift.matrix[h];
	  Ytrans = 1;
#ifdef DPD_DEBUG	  
	  yrow = Yt.shift.coltot[h];
	  ycol = Yt.shift.rowtot[h];
#endif	
	}
      else if(target_Y == 3) {
	  dpd_buf4_mat_irrep_init(Y, h);
	  dpd_buf4_mat_irrep_rd(Y, h);
	  dpd_buf4_mat_irrep_shift31(Y, h);
	  rking = 1;
	  Ymat = Y->shift.matrix[h];
	  Ytrans = 0;
#ifdef DPD_DEBUG	  
	  yrow = Y->shift.rowtot[h];
	  ycol = Y->shift.coltot[h];
#endif	  
	}
      else {
	  fprintf(stderr, "Junk Y index %d in contract442\n", target_Y);
	  exit(target_Y);
	}

      if(rking)
	  for(Gtar=0; Gtar < nirreps; Gtar++) {
#ifdef DPD_DEBUG
	      if((xrow[Gtar] != zrow[Gtar]) || (ycol[Gtar] != zcol[Gtar]) ||
		 (xcol[Gtar] != yrow[Gtar])) {
		  fprintf(stderr, "** Alignment error in contract442 **\n");
		  fprintf(stderr, "** Irrep %d; Subirrep %d **\n", h,Gtar);
		  dpd_error("dpd_contract442", stderr);
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
		  fprintf(stderr, "** Alignment error in contract442 **\n");
		  fprintf(stderr, "** Irrep %d; Subirrep %d **\n", h,Gtar);
		  dpd_error("dpd_contract442", stderr);
		}
#endif	      
	      newmm(Xmat[Gtar], Xtrans, Ymat[Gtar], Ytrans,
		    Z->matrix[Gtar], Z->params->rowtot[Gtar],
		    numlinks[Gtar], Z->params->coltot[Gtar],
		    alpha, 1.0);
	    }


      if(target_X == 0) dpd_buf4_mat_irrep_close(X, h);
      else if(target_X == 1) dpd_trans4_mat_irrep_close(&Xt, h);
      else if(target_X == 2) dpd_trans4_mat_irrep_close(&Xt, h);
      else if(target_X == 3) dpd_buf4_mat_irrep_close(X, h);

      if(target_Y == 0) dpd_buf4_mat_irrep_close(Y, h);
      else if(target_Y == 1) dpd_trans4_mat_irrep_close(&Yt, h);
      else if(target_Y == 2) dpd_trans4_mat_irrep_close(&Yt, h);
      else if(target_Y == 3) dpd_buf4_mat_irrep_close(Y, h);
    }

  if((target_X == 1) || (target_X == 2)) dpd_trans4_close(&Xt);
  if((target_Y == 1) || (target_Y == 2)) dpd_trans4_close(&Yt);


  dpd_file2_mat_wrt(Z);
  dpd_file2_mat_close(Z);

  return 0;
}


