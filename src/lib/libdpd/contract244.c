#include <stdio.h>
#include <math.h>
#include <libqt/qt.h>
#include "dpd.h"
#define EXTERN
#include "dpd.gbl"

/* dpd_contract244(): Contracts a two-index quantity with a
** four-index quantity to compute the contribution to another
** four-index quantity, beta * Z = alpha * X * Y.
**
** Arguments:
**   dpdfile2 *X: A pointer to the two-index file.
**   dpdbuf4 *Y: A pointer to the four-index buffer.
**   dpdbuf4 *Z: A pointer to the target four-index buffer.
**   int sum_X: Indicates which index on X is to be summed (takes a value of
**              0 or 1).
**   int sum_Y: Indicates which index on Y is to be summed (takes a value of
**              0, 1, 2, or 3).
**   int trans_Z: A boolean which indicates whether the indices in Z
**                must be transposed in order to match the external
**                index ordering in the product X * Y.
**   double alpha: A prefactor for the product alpha * X * Y.
**   double beta: A prefactor for the target beta * Z.
*/

int dpd_contract244(dpdfile2 *X, dpdbuf4 *Y, dpdbuf4 *Z, int sum_X, int sum_Y,
		    int trans_Z, double alpha, double beta)
{
  int h, nirreps, Gtar;
  int rking=0;
  int Xtrans,Ytrans;
  int *numlinks, *numrows, *numcols;
  dpdtrans4 Yt, Zt;
  double ***Ymat, ***Zmat;
#ifdef DPD_DEBUG
  int *xrow, *xcol, *yrow, *ycol, *zrow, *zcol;
#endif  

  nirreps = X->params->nirreps;

  dpd_file2_mat_init(X);
  dpd_file2_mat_rd(X);
  
  if(sum_X == 0) { Xtrans = 1; numlinks = X->params->rowtot; }
  else if(sum_X == 1) { Xtrans = 0; numlinks = X->params->coltot; }
  else { fprintf(stderr, "Junk X index %d\n", sum_X); exit(sum_X); }

  if((sum_Y == 1) || (sum_Y == 2)) dpd_trans4_init(&Yt, Y);

  if(trans_Z) dpd_trans4_init(&Zt, Z);

/*  if(fabs(beta) > 0.0) dpd_buf4_scm(Z, beta); */
  dpd_buf4_scm(Z, beta);

#ifdef DPD_DEBUG  
  if(Xtrans) { xrow = X->params->coltot; xcol = X->params->rowtot; }
  else { xrow = X->params->rowtot; xcol = X->params->coltot; }
#endif  

  for(h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(Z, h);
      if(fabs(beta) > 0.0) dpd_buf4_mat_irrep_rd(Z, h);
      if(trans_Z) {
	  dpd_trans4_mat_irrep_init(&Zt, h);
	  dpd_trans4_mat_irrep_rd(&Zt, h);
	  dpd_buf4_mat_irrep_close(Z, h);
	  dpd_trans4_mat_irrep_shift13(&Zt, h);
	  numrows = Zt.shift.rowtot[h];
	  numcols = Zt.shift.coltot[h];
	  Zmat = Zt.shift.matrix[h];
#ifdef DPD_DEBUG
	  zrow = Zt.shift.rowtot[h];
	  zcol = Zt.shift.coltot[h];
#endif
	}
      else {
	  dpd_buf4_mat_irrep_shift13(Z, h);
	  numrows = Z->shift.rowtot[h];
	  numcols = Z->shift.coltot[h];
	  Zmat = Z->shift.matrix[h];
#ifdef DPD_DEBUG
	  zrow = Z->shift.rowtot[h];
	  zcol = Z->shift.coltot[h];
#endif	  
	}

      if(sum_Y == 0) {
	  dpd_buf4_mat_irrep_init(Y, h);
	  dpd_buf4_mat_irrep_rd(Y, h);
	  dpd_buf4_mat_irrep_shift13(Y, h);
	  Ymat = Y->shift.matrix[h];
	  Ytrans = 0;
#ifdef DPD_DEBUG
	  yrow = Y->shift.rowtot[h];
	  ycol = Y->shift.coltot[h];
#endif
	}
      else if(sum_Y == 1) {
          dpd_buf4_mat_irrep_init(Y, h);
          dpd_buf4_mat_irrep_rd(Y, h);
          dpd_trans4_mat_irrep_init(&Yt, h);
          dpd_trans4_mat_irrep_rd(&Yt, h);
          dpd_buf4_mat_irrep_close(Y, h);
          dpd_trans4_mat_irrep_shift31(&Yt, h);
	  rking = 1;
          Ymat = Yt.shift.matrix[h];
          Ytrans = 1;
#ifdef DPD_DEBUG
	  yrow = Yt.shift.coltot[h];
	  ycol = Yt.shift.rowtot[h];
#endif  
        }
      else if(sum_Y == 2) {
          dpd_buf4_mat_irrep_init(Y, h);
          dpd_buf4_mat_irrep_rd(Y, h);
          dpd_trans4_mat_irrep_init(&Yt, h);
          dpd_trans4_mat_irrep_rd(&Yt, h);
          dpd_buf4_mat_irrep_close(Y, h);
          dpd_trans4_mat_irrep_shift13(&Yt, h);
          Ymat = Yt.shift.matrix[h];
          Ytrans = 0;
#ifdef DPD_DEBUG
	  yrow = Yt.shift.rowtot[h];
	  ycol = Yt.shift.coltot[h];
#endif 	  
        }
      else if(sum_Y == 3) {
          dpd_buf4_mat_irrep_init(Y, h);
          dpd_buf4_mat_irrep_rd(Y, h);
          dpd_buf4_mat_irrep_shift31(Y, h);
	  rking = 1;
          Ymat = Y->shift.matrix[h];
          Ytrans = 1;
#ifdef DPD_DEBUG
	  yrow = Y->shift.coltot[h];
	  ycol = Y->shift.rowtot[h];
#endif		  
        }

      if(rking) 
	  for(Gtar=0; Gtar < nirreps; Gtar++) {
#ifdef DPD_DEBUG
	      if((xrow[Gtar] != zrow[Gtar]) || (ycol[Gtar] != zcol[Gtar]) ||
		 (xcol[Gtar] != yrow[Gtar])) {
		  fprintf(stderr, "** Alignment error in contract244 **\n");
		  fprintf(stderr, "** Irrep: %d; Subirrep: %d **\n",h,Gtar);
		  dpd_error("dpd_contract244", stderr);
		}
#endif
	      newmm_rking(X->matrix[Gtar],Xtrans, Ymat[Gtar], Ytrans,
			  Zmat[Gtar], numrows[Gtar], numlinks[Gtar],
			  numcols[Gtar], alpha, 1.0);
	    }
      else
	  for(Gtar=0; Gtar < nirreps; Gtar++) {
#ifdef DPD_DEBUG
	      if((xrow[Gtar] != zrow[Gtar]) || (ycol[Gtar] != zcol[Gtar]) ||
		 (xcol[Gtar] != yrow[Gtar])) {
		  fprintf(stderr, "** Alignment error in contract244 **\n");
		  fprintf(stderr, "** Irrep: %d; Subirrep: %d **\n",h,Gtar);
		  dpd_error("dpd_contract244", stderr);
		}
#endif	      
	      newmm(X->matrix[Gtar],Xtrans, Ymat[Gtar], Ytrans,
		    Zmat[Gtar], numrows[Gtar], numlinks[Gtar],
		    numcols[Gtar], alpha, 1.0);
	    }

      if(sum_Y == 0) dpd_buf4_mat_irrep_close(Y, h);
      else if(sum_Y == 1) dpd_trans4_mat_irrep_close(&Yt, h);
      else if(sum_Y == 2) dpd_trans4_mat_irrep_close(&Yt, h);
      else if(sum_Y == 3) dpd_buf4_mat_irrep_close(Y, h);

      if(trans_Z) {
	  dpd_buf4_mat_irrep_init(Z, h);
	  dpd_trans4_mat_irrep_wrt(&Zt, h);
	  dpd_trans4_mat_irrep_close(&Zt, h);
	}

      dpd_buf4_mat_irrep_wrt(Z, h);
      dpd_buf4_mat_irrep_close(Z, h);

    }

  if((sum_Y == 1) || (sum_Y == 2)) dpd_trans4_close(&Yt);

  if(trans_Z) dpd_trans4_close(&Zt);

  dpd_file2_mat_close(X);

  return 0;
}
