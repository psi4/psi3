#include <stdio.h>
#include <math.h>
#include <qt.h>
#include "dpd.h"
#define EXTERN
#include "dpd.gbl"

/* dpd_contract221(): Contracts two-electron and one-electron quantities to
** give a product two-electron quantity.
**
** Arguments:
**   struct dpdbuf *X: A pointer to the leftmost dpd two-electron
**                     buffer in the product.
**   struct oe_dpdfile *Y: A pointer to the rightmost dpd one-electron
**                         file in the product.
**   struct dpdbuf *Z: A pointer to the dpd two-electron buffer target.
**   int sum_X: Indicates which index (values of 0, 1, 2, and 3) of X
**              is to be summed.
**   int sum_Y: Indicates which index (values of 0 and 1) of Y is to be summed.
**   int trans_Z: Boolean to indicate whether the final product must
**                be bra-ket transposed in order for its indices to
**                match those of the target, Z.
**   double alpha: A prefactor for the product alpha * X * Y.
**   double beta: A prefactor for the target beta * Z.
**   int print_flag: A boolean for the print routines.
**   FILE *outfile: The formatted output file stream.
*/

int dpd_contract221(struct dpdbuf *X, struct oe_dpdfile *Y, struct
		    dpdbuf *Z, int sum_X, int sum_Y, int trans_Z, double
		    alpha, double beta, int print_flag, FILE *outfile)
{
  int h, nirreps, Gtar;
  int rking=0;
  int Xtrans,Ytrans;
  int *numlinks, *numrows, *numcols;
  int incore, core, memoryd;
  int xcount, zcount, scount, Ysym;
  int rowx, rowz, colx, colz;
  int pq, rs, r, s, Gr, Gs;
  struct dpdtrans Xt, Zt;
  double ***Xmat, ***Zmat;
#ifdef DPD_DEBUG
  int *xrow, *xcol, *yrow, *ycol, *zrow, *zcol;
#endif  

  nirreps = X->params->nirreps;

  memoryd = dpd_memory/sizeof(double);
  incore = 1; /* default */

  dpd_oe_file_mat_init(Y);
  dpd_oe_file_mat_rd(Y, print_flag, outfile);
  if(sum_Y == 0) { Ytrans = 0; numlinks = Y->params->rowtot; }
  else if(sum_Y == 1) { Ytrans = 1; numlinks = Y->params->coltot; }
  else { fprintf(outfile, "Junk Y index %d\n", sum_Y); exit(sum_Y); }

  if((sum_X == 1) || (sum_X == 2)) dpd_trans_init(&Xt, X);

  if(trans_Z) dpd_trans_init(&Zt, Z);

  if(fabs(beta) > 0.0) dpd_scm(Z, beta, print_flag, outfile);

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
	  
      dpd_buf_mat_irrep_init(Z, h);
      if(fabs(beta) > 0.0) dpd_buf_mat_irrep_rd(Z, h, print_flag, outfile);
      if(trans_Z) {
	  dpd_trans_mat_irrep_init(&Zt, h);
	  dpd_trans_mat_irrep_rd(&Zt, h);
	  dpd_buf_mat_irrep_close(Z, h);
	  dpd_trans_mat_irrep_shift31(&Zt, h, print_flag, outfile);
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
	  dpd_buf_mat_irrep_shift31(Z, h, print_flag, outfile);
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
          dpd_buf_mat_irrep_init(X, h);
          dpd_buf_mat_irrep_rd(X, h, print_flag, outfile);
          dpd_buf_mat_irrep_shift13(X, h, print_flag, outfile);
          Xmat = X->shift.matrix[h];
          Xtrans = 1;
#ifdef DPD_DEBUG
	  xrow = X->shift.coltot[h];
	  xcol = X->shift.rowtot[h];
#endif	  
        }
      else if(sum_X == 1) {
          dpd_buf_mat_irrep_init(X, h);
          dpd_buf_mat_irrep_rd(X, h, print_flag, outfile);
          dpd_trans_mat_irrep_init(&Xt, h);
          dpd_trans_mat_irrep_rd(&Xt, h);
          dpd_buf_mat_irrep_close(X, h);
          dpd_trans_mat_irrep_shift31(&Xt, h, print_flag, outfile);
	  rking = 1;
          Xmat = Xt.shift.matrix[h];
          Xtrans = 0;
#ifdef DPD_DEBUG
	  xrow = Xt.shift.rowtot[h];
	  xcol = Xt.shift.coltot[h];
#endif 	  
        }
      else if(sum_X == 2) {
          dpd_buf_mat_irrep_init(X, h);
          dpd_buf_mat_irrep_rd(X, h, print_flag, outfile);
          dpd_trans_mat_irrep_init(&Xt, h);
          dpd_trans_mat_irrep_rd(&Xt, h);
          dpd_buf_mat_irrep_close(X, h);
          dpd_trans_mat_irrep_shift13(&Xt, h, print_flag, outfile);
          Xmat = Xt.shift.matrix[h];
          Xtrans = 1;
#ifdef DPD_DEBUG
	  xrow = Xt.shift.coltot[h];
	  xcol = Xt.shift.rowtot[h];
#endif	  
        }
      else if(sum_X == 3) {
          dpd_buf_mat_irrep_init(X, h);
          dpd_buf_mat_irrep_rd(X, h, print_flag, outfile);
          dpd_buf_mat_irrep_shift31(X, h, print_flag, outfile);
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
		  fprintf(outfile, "** Alignment error in contract221 **\n");
		  fprintf(outfile, "** Irrep: %d; Subirrep: %d **\n",h,Gtar);
		  dpd_error("dpd_contract221", outfile);
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
		  fprintf(outfile, "** Alignment error in contract221 **\n");
		  fprintf(outfile, "** Irrep: %d; Subirrep: %d **\n",h,Gtar);
		  dpd_error("dpd_contract221", outfile);
		}
#endif	    	      
	      newmm(Xmat[Gtar], Xtrans, Y->matrix[Gtar], Ytrans,
		    Zmat[Gtar], numrows[Gtar],numlinks[Gtar],
		    numcols[Gtar], alpha, 1.0);
	    }
      
      if(sum_X == 0) dpd_buf_mat_irrep_close(X, h);
      else if(sum_X == 1) dpd_trans_mat_irrep_close(&Xt, h);
      else if(sum_X == 2) dpd_trans_mat_irrep_close(&Xt, h);
      else if(sum_X == 3) dpd_buf_mat_irrep_close(X, h);

      if(trans_Z) {
	  dpd_buf_mat_irrep_init(Z, h);
	  dpd_trans_mat_irrep_wrt(&Zt, h);
	  dpd_trans_mat_irrep_close(&Zt, h);
	}

      dpd_buf_mat_irrep_wrt(Z, h, print_flag, outfile);
      dpd_buf_mat_irrep_close(Z, h);

    }  /* end if(incore) */
    else { /* out-of-core for "normal" 221 contractions */
	/* Prepare the input buffer for the X factor and the target*/
#ifdef DPD_DEBUG	
	fprintf(outfile, "\t221 out-of-core: %d\n", h);
#endif
	dpd_buf_mat_irrep_row_init(X, h);
	dpd_buf_mat_irrep_row_init(Z, h);
	
	/* Loop over rows of the X factor and the target */
	for(pq=0; pq < Z->params->rowtot[h]; pq++) {

	    dpd_buf_mat_irrep_row_zero(X, h, pq);
	    dpd_buf_mat_irrep_row_rd(X, h, pq, 0, outfile);

	    dpd_buf_mat_irrep_row_zero(Z, h, pq);

	    if(fabs(beta) > 0.0)
		dpd_buf_mat_irrep_row_rd(Z, h, pq, 0, outfile);

	    /*
	    for(rs=0; rs < X->params->coltot[h]; rs++)
		fprintf(outfile, "\t%d X[%d] = %20.15f\n", pq, rs,
			X->matrix[h][0][rs]);
			

	    for(rs=0; rs < Z->params->coltot[h]; rs++)
		fprintf(outfile, "\t%d Z[%d] = %20.15f\n", pq, rs,
			Z->matrix[h][0][rs]);
	    fflush(outfile);
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

	    dpd_buf_mat_irrep_row_wrt(Z, h, pq, 0, outfile);

	    /*
	    for(rs=0; rs < Z->params->coltot[h]; rs++)
		fprintf(outfile, "\t%d Zout[%d] = %20.15f\n", pq, rs,
			Z->matrix[h][0][rs]);
			*/

	  }

	dpd_buf_mat_irrep_row_close(X, h);
	dpd_buf_mat_irrep_row_close(Z, h);
      }
    }

  if((sum_X == 1) || (sum_X == 2)) dpd_trans_close(&Xt);

  if(trans_Z) dpd_trans_close(&Zt);

  dpd_oe_file_mat_close(Y);

  return 0;


}
