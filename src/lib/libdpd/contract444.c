#include <stdio.h>
#include <math.h>
#include <libqt/qt.h>
#include <libpsio/psio.h>
#include "dpd.h"
#define EXTERN
#include "dpd.gbl"

/* dpd_contract444(): Contracts a pair of four-index quantities to
** give a product four-index quantity.
**
** Arguments:
**   dpdbuf4 *X: A pointer to the leftmost dpd four-index
**               buffer in the product.
**   dpdbuf4 *Y: A pointer to the rightmost dpd four-index
**               buffer in the product.
**   int target_X: Indicates which pair of indices (0 = bra, 1 =
**                 ket) of X is the target pair.
**   int target_Y: Indicates which pair of indices (0 = bra, 1 =
**                 ket) of Y is the target pair.
**   double alpha: A prefactor for the product alpha * X * Y.
**   double beta: A prefactor for the target beta * Z.
*/

int dpd_contract444(dpdbuf4 *X, dpdbuf4 *Y, dpdbuf4 *Z,
		    int target_X, int target_Y, double alpha,
		    double beta)
{
  int n, hx, hy, hz, GX, GY, GZ, nirreps, Xtrans, Ytrans, *numlinks, symlink;
  int incore, memoryd, core, rows_per_bucket, nbuckets, rows_left, memtotal;
#ifdef DPD_DEBUG
  int *xrow, *xcol, *yrow, *ycol, *zrow, *zcol;
  double byte_conv;
#endif

  nirreps = X->params->nirreps;
  GX = X->file.my_irrep;
  GY = Y->file.my_irrep;
  GZ = Z->file.my_irrep;

  if(target_X == 0) { Xtrans = 0; numlinks = X->params->coltot; symlink=GX; }
  else if(target_X == 1) { Xtrans = 1; numlinks = X->params->rowtot; symlink=0; }

  if(target_Y == 0) Ytrans = 1;
  else if(target_Y == 1) Ytrans = 0;

#ifdef DPD_DEBUG
  if(Xtrans) { xrow = X->params->coltot; xcol = X->params->rowtot; }
  else { xrow = X->params->rowtot; xcol = X->params->coltot; }

  if(Ytrans) { yrow = Y->params->coltot; ycol = Y->params->rowtot; }
  else { yrow = Y->params->rowtot; ycol = Y->params->coltot; }

  zrow = Z->params->rowtot; zcol = Z->params->coltot;
  
  if((zrow != xrow) || (zcol != ycol) || (xcol != yrow)) {
    fprintf(stderr, "** Alignment error in contract444 **\n");
    dpd_error("dpd_contract444",stderr);
  }

#endif
  

  for(hx=0; hx < nirreps; hx++) {

    if      ((!Xtrans)&&(!Ytrans))  {hy = hx^GX;    hz = hx;    }
    else if ((!Xtrans)&&( Ytrans))  {hy = hx^GX^GY; hz = hx;    }
    else if (( Xtrans)&&(!Ytrans))  {hy = hx;       hz = hx^GX; }
    else /* (( Xtrans)&&( Ytrans))*/{hy = hx^GY;    hz = hx^GX; }

    dpd_buf4_mat_irrep_init(Y, hy);
    dpd_buf4_mat_irrep_rd(Y, hy);
    dpd_buf4_mat_irrep_init(Z, hz);
    if(fabs(beta) > 0.0) 
      dpd_buf4_mat_irrep_rd(Z, hz);
	
    memoryd = dpd_memfree();

    if(X->params->rowtot[hx] && X->params->coltot[hx^GX]) {

      if(X->params->coltot[hx^GX])
	rows_per_bucket = memoryd/X->params->coltot[hx^GX];
      else rows_per_bucket = -1;

      if(rows_per_bucket > X->params->rowtot[hx])
	rows_per_bucket = X->params->rowtot[hx];

      if(!rows_per_bucket)
	dpd_error("contract444: Not enough memory for one row", stderr);

      nbuckets = ceil((double) X->params->rowtot[hx]/
		      (double) rows_per_bucket);

      rows_left = X->params->rowtot[hx] % rows_per_bucket;
      
      incore = 1;
      if(nbuckets > 1) incore = 0;
    }
    else incore = 1;

#ifdef DPD_DEBUG
    fprintf(stderr, "Contract444: memory information.\n");
    fprintf(stderr, "Contract444: h = %d, row = %d, col = %d, tot = %d\n", 
            h, X->params->rowtot[hx], X->params->coltot[hx^GX],
            X->params->rowtot[hx] * X->params->coltot[hx^GX]);
    if(!incore) {
      fprintf(stderr, "Contract444: nbuckets = %d\n", nbuckets);
      fprintf(stderr, "Contract444: rows_per_bucket = %d\n",rows_per_bucket);
      fprintf(stderr, "Contract444: rows_left = %d\n",rows_left);
      memtotal = X->params->rowtot[hx] * X->params->coltot[hx^GX];
      byte_conv = ((double) sizeof(double))/1e6;
      fprintf(stderr, "Contract444: out of core algorithm used.\n");
      fprintf(stderr, "Contract444: memtotal = %d.\n", memtotal);
      fprintf(stderr, "Contract444: Need %5.2f MB to run in memory.\n",
	      ((double) memtotal)*byte_conv);
      dpd_file4_cache_print(stderr);
    }
    fflush(stderr);
#endif

    if(!incore && Xtrans) {
      dpd_file4_cache_print(stderr);
      dpd_error("out-of-core contract444 Xtrans=1 not coded", stderr);
    }

    if(incore) {
      if(fabs(beta) > 0.0) dpd_buf4_mat_irrep_rd(Z, hz);
	  
      dpd_buf4_mat_irrep_init(X, hx);
      dpd_buf4_mat_irrep_rd(X, hx);

      newmm(X->matrix[hx], Xtrans, Y->matrix[hy], Ytrans,
	    Z->matrix[hz], Z->params->rowtot[hz], numlinks[hx^symlink],
	    Z->params->coltot[hz^GZ], alpha, beta);

      dpd_buf4_mat_irrep_close(X, hx);

      dpd_buf4_mat_irrep_wrt(Z, hz);

    }
    else {

      if(!Ytrans) { 
	fprintf(stderr, "Contract444: Problem!  Out-of-core algorithm used, but Ytrans == 0!\n");
	exit(2);
      }

      dpd_buf4_mat_irrep_init_block(X, hx, rows_per_bucket);

      for(n=0; n < (rows_left ? nbuckets-1 : nbuckets); n++) {

	dpd_buf4_mat_irrep_rd_block(X, hx, n*rows_per_bucket, rows_per_bucket);

	C_DGEMM('n', 't', rows_per_bucket, Z->params->coltot[hz^GZ],
		numlinks[hx^symlink], alpha, &(X->matrix[hx][0][0]), numlinks[hx^symlink],
		&(Y->matrix[hy][0][0]), numlinks[hx^symlink], beta,
		&(Z->matrix[hz][n*rows_per_bucket][0]), Z->params->coltot[hz^GZ]);

      }

      if(rows_left) {

	dpd_buf4_mat_irrep_rd_block(X, hx, n*rows_per_bucket, rows_left);

	C_DGEMM('n', 't', rows_left, Z->params->coltot[hz^GZ],
		numlinks[hx^symlink], alpha, &(X->matrix[hx][0][0]), numlinks[hx^symlink],
		&(Y->matrix[hy][0][0]), numlinks[hx^symlink], beta,
		&(Z->matrix[hz][n*rows_per_bucket][0]), Z->params->coltot[hz^GZ]);

      }
	      
      dpd_buf4_mat_irrep_close_block(X, hx, rows_per_bucket);

    }

    dpd_buf4_mat_irrep_close(Y, hy);
    dpd_buf4_mat_irrep_wrt(Z, hz);
    dpd_buf4_mat_irrep_close(Z, hz);
  }

  return 0;
}
