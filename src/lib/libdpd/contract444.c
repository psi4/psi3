#include <stdio.h>
#include <math.h>
#include <qt.h>
#include <psio.h>
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
**   int target_X: Indicates the which pair of indices (0 = bra, 1 =
**                 ket) of X is the target pair.
**   int target_Y: Indicates the which pair of indices (0 = bra, 1 =
**                 ket) of Y is the target pair.
**   double alpha: A prefactor for the product alpha * X * Y.
**   double beta: A prefactor for the target beta * Z.
*/

int dpd_contract444(dpdbuf4 *X, dpdbuf4 *Y, dpdbuf4 *Z,
		    int target_X, int target_Y, double alpha,
		    double beta)
{
  int n, h, nirreps, Xtrans, Ytrans, *numlinks;
  int incore, memoryd, core, rows_per_bucket, nbuckets, rows_left, memtotal;
#ifdef DPD_DEBUG
  int *xrow, *xcol, *yrow, *ycol, *zrow, *zcol;
#endif

  nirreps = X->params->nirreps;

  if(target_X == 0) { Xtrans = 0; numlinks = X->params->coltot; }
  else if(target_X == 1) { Xtrans = 1; numlinks = X->params->rowtot; }

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
  
  memoryd = dpd_default->memory;

  for(h=0; h < nirreps; h++) {

      /* Compute the memory requirements for the Z and Y */
      core = Z->params->rowtot[h] * Z->params->coltot[h] +
	     Y->params->rowtot[h] * Y->params->coltot[h];

      if(!X->params->rowtot[h]) rows_per_bucket = 0;
      else rows_per_bucket = (memoryd-core)/X->params->rowtot[h];
      if(rows_per_bucket > X->params->rowtot[h])
	  rows_per_bucket = X->params->rowtot[h];

      if(!rows_per_bucket) nbuckets = 1;
      else
        nbuckets = 
           ceil((double) X->params->rowtot[h]/(double) rows_per_bucket);
      if(!rows_per_bucket) rows_left = 0;
      else rows_left = X->params->rowtot[h] % rows_per_bucket;
      
      incore = 1;
      if(nbuckets > 1) incore = 0;

#ifdef DPD_DEBUG
      if(!incore) {
	  memtotal = core + X->params->rowtot[h] * X->params->coltot[h];
	  memtotal *= sizeof(double);
	  fprintf(stderr, "Contract444: out of core algorithm used.\n");
	  fprintf(stderr, "Need %5.2f MB to run in memory.\n",
		  ((double) memtotal)/1e6);
	}
#endif

      if(!incore && Xtrans) {
	dpd_file4_cache_print(stderr);
	dpd_error("out-of-core contract444 Xtrans=1 not coded", stderr);
	}

      dpd_buf4_mat_irrep_init(Y, h);
      dpd_buf4_mat_irrep_rd(Y, h);
	
      if(incore) {

	  dpd_buf4_mat_irrep_init(Z, h);
	  
	  if(fabs(beta) > 0.0) dpd_buf4_mat_irrep_rd(Z, h);
	  
	  dpd_buf4_mat_irrep_init(X, h);
	  dpd_buf4_mat_irrep_rd(X, h);

	  newmm(X->matrix[h], Xtrans, Y->matrix[h], Ytrans,
		Z->matrix[h], Z->params->rowtot[h], numlinks[h],
		Z->params->coltot[h], alpha, beta);

	  dpd_buf4_mat_irrep_close(X, h);

	  dpd_buf4_mat_irrep_wrt(Z, h);
	  dpd_buf4_mat_irrep_close(Z, h);

	}
      else {

	  dpd_buf4_mat_irrep_init_block(X, h, rows_per_bucket);
	  dpd_buf4_mat_irrep_init_block(Z, h, rows_per_bucket);

	  for(n=0; n < (rows_left ? nbuckets-1 : nbuckets); n++) {
	      
	      dpd_buf4_mat_irrep_rd_block(X, h, n*rows_per_bucket,
					 rows_per_bucket);

	      if(fabs(beta) > 0.0)
		  dpd_buf4_mat_irrep_rd_block(Z, h, n*rows_per_bucket,
					     rows_per_bucket);
	      
	      C_DGEMM('n', 't', rows_per_bucket, Z->params->coltot[h],
		      numlinks[h], alpha, &(X->matrix[h][0][0]), numlinks[h],
		      &(Y->matrix[h][0][0]), numlinks[h], beta,
		      &(Z->matrix[h][0][0]), Z->params->coltot[h]);

	      dpd_buf4_mat_irrep_wrt_block(Z, h, n*rows_per_bucket,
					  rows_per_bucket);
	    }

	if(rows_left) {

	    dpd_buf4_mat_irrep_rd_block(X, h, n*rows_per_bucket, rows_left);

	    if(fabs(beta) > 0.0)
		dpd_buf4_mat_irrep_rd_block(Z, h, n*rows_per_bucket, rows_left);
	      
	    C_DGEMM('n', 't', rows_left, Z->params->coltot[h],
		    numlinks[h], alpha, &(X->matrix[h][0][0]), numlinks[h],
		    &(Y->matrix[h][0][0]), numlinks[h], beta,
		    &(Z->matrix[h][0][0]), Z->params->coltot[h]);

	    dpd_buf4_mat_irrep_wrt_block(Z, h, n*rows_per_bucket, rows_left);
	  }
	      
	dpd_buf4_mat_irrep_close_block(X, h, rows_per_bucket);
	dpd_buf4_mat_irrep_close_block(Z, h, rows_per_bucket);

	}

      dpd_buf4_mat_irrep_close(Y, h);
    }

  return 0;
}
