#include <stdio.h>
#include <math.h>
#include <libqt/qt.h>
#include "dpd.h"
#define EXTERN
#include "dpd.gbl"

/* dpd_contract424(): Contracts four-index and two-index quantities to
** give a product four-index quantity.
**
** Arguments:
**   dpdbuf4 *X: A pointer to the leftmost dpd four-index
**               buffer in the product.
**   dpdfile2 *Y: A pointer to the rightmost dpd two-index
**                file in the product.
**   dpdbuf4 *Z: A pointer to the dpd four-index buffer target.
**   int sum_X: Indicates which index (values of 0, 1, 2, and 3) of X
**              is to be summed.
**   int sum_Y: Indicates which index (values of 0 and 1) of Y is to be summed.
**   int trans_Z: Boolean to indicate whether the final product must
**                be bra-ket transposed in order for its indices to
**                match those of the target, Z.
**   double alpha: A prefactor for the product alpha * X * Y.
**   double beta: A prefactor for the target beta * Z.
*/

int dpd_contract424(dpdbuf4 *X, dpdfile2 *Y, dpdbuf4 *Z, int sum_X,
		    int sum_Y, int trans_Z, double alpha, double beta)
{
  int h, nirreps, Gtar;
  int rking=0;
  int Xtrans,Ytrans;
  int *numlinks, *numrows, *numcols;
  int incore, core, memoryd;
  int xcount, zcount, scount, Ysym;
  int rowx, rowz, colx, colz;
  int pq, rs, r, s, Gr, Gs;
  dpdtrans4 Xt, Zt;
  double ***Xmat, ***Zmat;
#ifdef DPD_DEBUG
  int *xrow, *xcol, *yrow, *ycol, *zrow, *zcol;
#endif  

  nirreps = X->params->nirreps;

  memoryd = dpd_main.memory;
  incore = 1; /* default */

  dpd_file2_mat_init(Y);
  dpd_file2_mat_rd(Y);
  if(sum_Y == 0) { Ytrans = 0; numlinks = Y->params->rowtot; }
  else if(sum_Y == 1) { Ytrans = 1; numlinks = Y->params->coltot; }
  else { fprintf(stderr, "Junk Y index %d\n", sum_Y); exit(sum_Y); }

  if((sum_X == 1) || (sum_X == 2)) dpd_trans4_init(&Xt, X);

  if(trans_Z) dpd_trans4_init(&Zt, Z);

/*  if(fabs(beta) > 0.0) dpd_buf4_scm(Z, beta); */
  dpd_buf4_scm(Z, beta);

#ifdef DPD_DEBUG  
  if(Ytrans) { yrow = Y->params->coltot; ycol = Y->params->rowtot; }
  else { yrow = Y->params->rowtot; ycol = Y->params->coltot; }
#endif  

  for(h=0; h < nirreps; h++) {

      /* Compute the core requirements for the straight contraction */
      core = Z->params->rowtot[h] * Z->params->coltot[h] +
	     X->params->rowtot[h] * X->params->coltot[h];

      /* Force incore for all but a "normal" 221 contraction for now */
      if(core > memoryd) incore = 0;
      if(trans_Z || sum_X == 0 || sum_X == 1 || sum_X == 2) incore = 1;

      if(incore) { 
	  
      dpd_buf4_mat_irrep_init(Z, h);
      if(fabs(beta) > 0.0) dpd_buf4_mat_irrep_rd(Z, h);
      if(trans_Z) {
	  dpd_trans4_mat_irrep_init(&Zt, h);
	  dpd_trans4_mat_irrep_rd(&Zt, h);
	  dpd_buf4_mat_irrep_close(Z, h);
	  dpd_trans4_mat_irrep_shift31(&Zt, h);
	  rking = 1;
	  numrows = Zt.shift.rowtot[h];
	  numcols = Zt.shift.coltot[h];
	  Zmat = Zt.shift.matrix[h];
#ifdef DPD_DEBUG
	  zrow = Zt.shift.rowtot[h];
	  zcol = Zt.shift.coltot[h];
#endif	  
	}
      else {
	  dpd_buf4_mat_irrep_shift31(Z, h);
	  rking = 1;
	  numrows = Z->shift.rowtot[h];
	  numcols = Z->shift.coltot[h];
	  Zmat = Z->shift.matrix[h];
#ifdef DPD_DEBUG
	  zrow = Z->shift.rowtot[h];
	  zcol = Z->shift.coltot[h];
#endif	  
	}

      if(sum_X == 0) {
          dpd_buf4_mat_irrep_init(X, h);
          dpd_buf4_mat_irrep_rd(X, h);
          dpd_buf4_mat_irrep_shift13(X, h);
          Xmat = X->shift.matrix[h];
          Xtrans = 1;
#ifdef DPD_DEBUG
	  xrow = X->shift.coltot[h];
	  xcol = X->shift.rowtot[h];
#endif	  
        }
      else if(sum_X == 1) {
          dpd_buf4_mat_irrep_init(X, h);
          dpd_buf4_mat_irrep_rd(X, h);
          dpd_trans4_mat_irrep_init(&Xt, h);
          dpd_trans4_mat_irrep_rd(&Xt, h);
          dpd_buf4_mat_irrep_close(X, h);
          dpd_trans4_mat_irrep_shift31(&Xt, h);
	  rking = 1;
          Xmat = Xt.shift.matrix[h];
          Xtrans = 0;
#ifdef DPD_DEBUG
	  xrow = Xt.shift.rowtot[h];
	  xcol = Xt.shift.coltot[h];
#endif 	  
        }
      else if(sum_X == 2) {
          dpd_buf4_mat_irrep_init(X, h);
          dpd_buf4_mat_irrep_rd(X, h);
          dpd_trans4_mat_irrep_init(&Xt, h);
          dpd_trans4_mat_irrep_rd(&Xt, h);
          dpd_buf4_mat_irrep_close(X, h);
          dpd_trans4_mat_irrep_shift13(&Xt, h);
          Xmat = Xt.shift.matrix[h];
          Xtrans = 1;
#ifdef DPD_DEBUG
	  xrow = Xt.shift.coltot[h];
	  xcol = Xt.shift.rowtot[h];
#endif	  
        }
      else if(sum_X == 3) {
          dpd_buf4_mat_irrep_init(X, h);
          dpd_buf4_mat_irrep_rd(X, h);
          dpd_buf4_mat_irrep_shift31(X, h);
	  rking = 1;
          Xmat = X->shift.matrix[h];
          Xtrans = 0;
#ifdef DPD_DEBUG
	  xrow = X->shift.rowtot[h];
	  xcol = X->shift.coltot[h];
#endif  
        }

      if(rking)
	  for(Gtar=0; Gtar < nirreps; Gtar++) {
#ifdef DPD_DEBUG
	      if((xrow[Gtar] != zrow[Gtar]) || (ycol[Gtar] != zcol[Gtar]) ||
		 (xcol[Gtar] != yrow[Gtar])) {
		  fprintf(stderr, "** Alignment error in contract424 **\n");
		  fprintf(stderr, "** Irrep: %d; Subirrep: %d **\n",h,Gtar);
		  dpd_error("dpd_contract424", stderr);
		}
#endif      
	      newmm_rking(Xmat[Gtar], Xtrans, Y->matrix[Gtar], Ytrans,
			  Zmat[Gtar], numrows[Gtar],numlinks[Gtar],
			  numcols[Gtar], alpha, 1.0);
	    }
      else 
	  for(Gtar=0; Gtar < nirreps; Gtar++) {
#ifdef DPD_DEBUG
	      if((xrow[Gtar] != zrow[Gtar]) || (ycol[Gtar] != zcol[Gtar]) ||
		 (xcol[Gtar] != yrow[Gtar])) {
		  fprintf(stderr, "** Alignment error in contract424 **\n");
		  fprintf(stderr, "** Irrep: %d; Subirrep: %d **\n",h,Gtar);
		  dpd_error("dpd_contract424", stderr);
		}
#endif	    	      
	      newmm(Xmat[Gtar], Xtrans, Y->matrix[Gtar], Ytrans,
		    Zmat[Gtar], numrows[Gtar],numlinks[Gtar],
		    numcols[Gtar], alpha, 1.0);
	    }
      
      if(sum_X == 0) dpd_buf4_mat_irrep_close(X, h);
      else if(sum_X == 1) dpd_trans4_mat_irrep_close(&Xt, h);
      else if(sum_X == 2) dpd_trans4_mat_irrep_close(&Xt, h);
      else if(sum_X == 3) dpd_buf4_mat_irrep_close(X, h);

      if(trans_Z) {
	  dpd_buf4_mat_irrep_init(Z, h);
	  dpd_trans4_mat_irrep_wrt(&Zt, h);
	  dpd_trans4_mat_irrep_close(&Zt, h);
	}

      dpd_buf4_mat_irrep_wrt(Z, h);
      dpd_buf4_mat_irrep_close(Z, h);

    }  /* end if(incore) */
    else { /* out-of-core for "normal" 424 contractions */
	/* Prepare the input buffer for the X factor and the target*/
#ifdef DPD_DEBUG	
	fprintf(stderr, "\t424 out-of-core: %d\n", h);
#endif
	dpd_buf4_mat_irrep_row_init(X, h);
	dpd_buf4_mat_irrep_row_init(Z, h);
	
	/* Loop over rows of the X factor and the target */
	for(pq=0; pq < Z->params->rowtot[h]; pq++) {

	    dpd_buf4_mat_irrep_row_zero(X, h, pq);
	    dpd_buf4_mat_irrep_row_rd(X, h, pq);

	    dpd_buf4_mat_irrep_row_zero(Z, h, pq);

	    if(fabs(beta) > 0.0)
		dpd_buf4_mat_irrep_row_rd(Z, h, pq);

	    /*
	    for(rs=0; rs < X->params->coltot[h]; rs++)
		fprintf(stderr, "\t%d X[%d] = %20.15f\n", pq, rs,
			X->matrix[h][0][rs]);
			

	    for(rs=0; rs < Z->params->coltot[h]; rs++)
		fprintf(stderr, "\t%d Z[%d] = %20.15f\n", pq, rs,
			Z->matrix[h][0][rs]);
	    fflush(stderr);
	    */

	    xcount = zcount = 0;

	    for(Gr=0; Gr < nirreps; Gr++) {
		Gs = Gr^h;

		rowx = X->params->rpi[Gr];
		colx = X->params->spi[Gs];
		rowz = Z->params->rpi[Gr];
		colz = Z->params->spi[Gs];

		if(rowx && colx && colz)
		    newmm2(&(X->matrix[h][0][xcount]),0,
			   &(Y->matrix[Gs][0][0]),Ytrans,
			   &(Z->matrix[h][0][zcount]),
			   rowx,colx,colz,alpha,1.0);

		xcount += rowx * colx;
		zcount += rowz * colz;

	      }

	    dpd_buf4_mat_irrep_row_wrt(Z, h, pq);

	    /*
	    for(rs=0; rs < Z->params->coltot[h]; rs++)
		fprintf(stderr, "\t%d Zout[%d] = %20.15f\n", pq, rs,
			Z->matrix[h][0][rs]);
			*/

	  }

	dpd_buf4_mat_irrep_row_close(X, h);
	dpd_buf4_mat_irrep_row_close(Z, h);
      }
    }

  if((sum_X == 1) || (sum_X == 2)) dpd_trans4_close(&Xt);

  if(trans_Z) dpd_trans4_close(&Zt);

  dpd_file2_mat_close(Y);

  return 0;


}
