#include <stdio.h>
#include <math.h>
#include <libqt/qt.h>
#include "dpd.h"
#define EXTERN
#include "dpd.gbl"
extern FILE *outfile;

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
 **   int Ztrans: Boolean to indicate whether the final product must
 **                be bra-ket transposed in order for its indices to
 **                match those of the target, Z.
 **   double alpha: A prefactor for the product alpha * X * Y.
 **   double beta: A prefactor for the target beta * Z.
 */

int dpd_contract424(dpdbuf4 *X, dpdfile2 *Y, dpdbuf4 *Z, int sum_X,
    int sum_Y, int Ztrans, double alpha, double beta)
{
  int h, nirreps, GX, GY, GZ, hxbuf, hzbuf, h0, Hx, Hy, Hz;
  int rking=0, symlink;
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
  GX = X->file.my_irrep; 
  GY = Y->my_irrep;
  GZ = Z->file.my_irrep;

  memoryd = dpd_main.memory;
  incore = 1; /* default */

  dpd_file2_mat_init(Y);
  dpd_file2_mat_rd(Y);
  if(sum_Y == 0) { Ytrans = 0; numlinks = Y->params->rowtot; symlink=0; }
  else if(sum_Y == 1) { Ytrans = 1; numlinks = Y->params->coltot; symlink=GY; }
  else { fprintf(stderr, "Junk Y index %d\n", sum_Y); exit(sum_Y); }

  if((sum_X == 1) || (sum_X == 2)) dpd_trans4_init(&Xt, X);

  if(Ztrans) dpd_trans4_init(&Zt, Z);

  /*  if(fabs(beta) > 0.0) dpd_buf4_scm(Z, beta); */
  dpd_buf4_scm(Z, beta);

#ifdef DPD_DEBUG  
  if(Ytrans) { yrow = Y->params->coltot; ycol = Y->params->rowtot; }
  else { yrow = Y->params->rowtot; ycol = Y->params->coltot; }
#endif  

  for(hxbuf=0; hxbuf < nirreps; hxbuf++) {

    if (sum_X < 2) {
      if (!Ztrans) hzbuf = hxbuf^GX; else hzbuf = hxbuf^GX^GZ;
    }
    else{
      if (!Ztrans) hzbuf = hxbuf; else hzbuf = hxbuf^GZ;
    }

    /* Compute the core requirements for the straight contraction */
    core = Z->params->rowtot[hzbuf] * Z->params->coltot[hzbuf^GZ] +
      X->params->rowtot[hxbuf] * X->params->coltot[hxbuf^GX];

    /* Force incore for all but a "normal" 221 contraction for now */
    if(core > memoryd) incore = 0;
    if(Ztrans || sum_X == 0 || sum_X == 1 || sum_X == 2) incore = 1;

    if(incore) { 

      dpd_buf4_mat_irrep_init(Z, hzbuf);
      if(fabs(beta) > 0.0) dpd_buf4_mat_irrep_rd(Z, hzbuf);
      if(Ztrans) {
        dpd_trans4_mat_irrep_init(&Zt, hzbuf);
        dpd_trans4_mat_irrep_rd(&Zt, hzbuf);
        dpd_buf4_mat_irrep_close(Z, hzbuf);
        dpd_trans4_mat_irrep_shift31(&Zt, hzbuf);
        rking = 1;
        numrows = Zt.shift.rowtot[hzbuf];
        numcols = Zt.shift.coltot[hzbuf];
        Zmat = Zt.shift.matrix[hzbuf];
#ifdef DPD_DEBUG
        zrow = Zt.shift.rowtot[hzbuf];
        zcol = Zt.shift.coltot[hzbuf];
#endif	  
      }
      else {
        dpd_buf4_mat_irrep_shift31(Z, hzbuf);
        rking = 1;
        numrows = Z->shift.rowtot[hzbuf];
        numcols = Z->shift.coltot[hzbuf];
        Zmat = Z->shift.matrix[hzbuf];
#ifdef DPD_DEBUG
        zrow = Z->shift.rowtot[hzbuf];
        zcol = Z->shift.coltot[hzbuf];
#endif	  
      }

      if(sum_X == 0) {
        dpd_buf4_mat_irrep_init(X, hxbuf);
        dpd_buf4_mat_irrep_rd(X, hxbuf);
        dpd_buf4_mat_irrep_shift13(X, hxbuf);
        Xmat = X->shift.matrix[hxbuf];
        Xtrans = 1;
#ifdef DPD_DEBUG
        xrow = X->shift.coltot[hxbuf];
        xcol = X->shift.rowtot[hxbuf];
#endif	  
      }
      else if(sum_X == 1) {
        dpd_buf4_mat_irrep_init(X, hxbuf);
        dpd_buf4_mat_irrep_rd(X, hxbuf);
        dpd_trans4_mat_irrep_init(&Xt, hxbuf);
        dpd_trans4_mat_irrep_rd(&Xt, hxbuf);
        dpd_buf4_mat_irrep_close(X, hxbuf);
        dpd_trans4_mat_irrep_shift31(&Xt, hxbuf);
        rking = 1;
        Xmat = Xt.shift.matrix[hxbuf];
        Xtrans = 0;
#ifdef DPD_DEBUG
        xrow = Xt.shift.rowtot[hxbuf];
        xcol = Xt.shift.coltot[hxbuf];
#endif 	  
      }
      else if(sum_X == 2) {
        dpd_buf4_mat_irrep_init(X, hxbuf);
        dpd_buf4_mat_irrep_rd(X, hxbuf);
        dpd_trans4_mat_irrep_init(&Xt, hxbuf);
        dpd_trans4_mat_irrep_rd(&Xt, hxbuf);
        dpd_buf4_mat_irrep_close(X, hxbuf);
        dpd_trans4_mat_irrep_shift13(&Xt, hxbuf);
        Xmat = Xt.shift.matrix[hxbuf];
        Xtrans = 1;
#ifdef DPD_DEBUG
        xrow = Xt.shift.coltot[hxbuf];
        xcol = Xt.shift.rowtot[hxbuf];
#endif	  
      }
      else if(sum_X == 3) {
        dpd_buf4_mat_irrep_init(X, hxbuf);
        dpd_buf4_mat_irrep_rd(X, hxbuf);
        dpd_buf4_mat_irrep_shift31(X, hxbuf);
        rking = 1;
        Xmat = X->shift.matrix[hxbuf];
        Xtrans = 0;
#ifdef DPD_DEBUG
        xrow = X->shift.rowtot[hxbuf];
        xcol = X->shift.coltot[hxbuf];
#endif  
      }

      if(rking)
        for(Hz=0; Hz < nirreps; Hz++) {
          if      (!Xtrans && !Ytrans) {Hx=Hz;       Hy = Hz^GX; }
          else if (!Xtrans &&  Ytrans) {Hx=Hz;       Hy = Hz^GX^GY; }
          else if ( Xtrans && !Ytrans) {Hx=Hz^GX;    Hy = Hz^GX; }
          else if ( Xtrans &&  Ytrans) {Hx=Hz^GX;    Hy = Hz^GX^GY; }
#ifdef DPD_DEBUG
          if((xrow[Hz] != zrow[Hz]) || (ycol[Hz] != zcol[Hz]) ||
              (xcol[Hz] != yrow[Hz])) {
            fprintf(stderr, "** Alignment error in contract424 **\n");
            fprintf(stderr, "** Irrep: %d; Subirrep: %d **\n",hxbuf,Hz);
            dpd_error("dpd_contract424", stderr);
          }
#endif      
	  /* fprintf(outfile,"Hx %d Hy %d Hz %d\n",Hx,Hy,Hz);
	     fprintf(outfile,"numrows %d numlinks %d numcols %d\n",numrows[Hz],numlinks[Hy],numcols[Hz]); */
          newmm_rking(Xmat[Hx], Xtrans, Y->matrix[Hy], Ytrans,
              Zmat[Hz], numrows[Hz], numlinks[Hy^symlink],
              numcols[Hz], alpha, 1.0);
        }
      else 
        for(Hz=0; Hz < nirreps; Hz++) {
          if      (!Xtrans && !Ytrans) {Hx=Hz;       Hy = Hz^GX; }
          else if (!Xtrans &&  Ytrans) {Hx=Hz;       Hy = Hz^GX^GY; }
          else if ( Xtrans && !Ytrans) {Hx=Hz^GX;    Hy = Hz^GX; }
          else if ( Xtrans &&  Ytrans) {Hx=Hz^GX;    Hy = Hz^GX^GY; }
#ifdef DPD_DEBUG
          if((xrow[Hz] != zrow[Hz]) || (ycol[Hz] != zcol[Hz]) ||
              (xcol[Hz] != yrow[Hz])) {
            fprintf(stderr, "** Alignment error in contract424 **\n");
            fprintf(stderr, "** Irrep: %d; Subirrep: %d **\n",hxbuf,Hz);
            dpd_error("dpd_contract424", stderr);
          }
#endif	    	      
	  if(numrows[Hz] && numcols[Hz] && numlinks[Hy^symlink]) {
	  if(!Xtrans && !Ytrans) {
	    C_DGEMM('n','n',numrows[Hz], numcols[Hz], numlinks[Hy^symlink],
		alpha, &(Xmat[Hz][0][0]), numlinks[Hy^symlink],
		&(Y->matrix[Hy][0][0]), numcols[Hz], 1.0,
		&(Zmat[Hz][0][0]), numcols[Hz]);
	  }
	  else if(Xtrans && !Ytrans) {
	    C_DGEMM('t','n',numrows[Hz], numcols[Hz], numlinks[Hy^symlink],
		alpha, &(Xmat[Hz][0][0]), numrows[Hz],
		&(Y->matrix[Hy][0][0]), numcols[Hz], 1.0,
		&(Zmat[Hz][0][0]), numcols[Hz]);
	  }
	  else if(!Xtrans && Ytrans) {
	    C_DGEMM('n','t',numrows[Hz], numcols[Hz], numlinks[Hy^symlink],
		alpha, &(Xmat[Hz][0][0]), numlinks[Hy^symlink],
		&(Y->matrix[Hy][0][0]), numlinks[Hy^symlink], 1.0,
		&(Zmat[Hz][0][0]), numcols[Hz]);
	  }
	  else {
	    C_DGEMM('t','t',numrows[Hz], numcols[Hz], numlinks[Hy^symlink],
		alpha, &(Xmat[Hz][0][0]), numrows[Hz],
		&(Y->matrix[Hy][0][0]), numlinks[Hy^symlink], 1.0,
		&(Zmat[Hz][0][0]), numcols[Hz]);
	  }
	  }
	  /*
          newmm(Xmat[Hx], Xtrans, Y->matrix[Hy], Ytrans,
              Zmat[Hz], numrows[Hz], numlinks[Hy^symlink],
              numcols[Hz], alpha, 1.0);
	      */
        }

      if(sum_X == 0) dpd_buf4_mat_irrep_close(X, hxbuf);
      else if(sum_X == 1) dpd_trans4_mat_irrep_close(&Xt, hxbuf);
      else if(sum_X == 2) dpd_trans4_mat_irrep_close(&Xt, hxbuf);
      else if(sum_X == 3) dpd_buf4_mat_irrep_close(X, hxbuf);

      if(Ztrans) {
        dpd_buf4_mat_irrep_init(Z, hzbuf);
        dpd_trans4_mat_irrep_wrt(&Zt, hzbuf);
        dpd_trans4_mat_irrep_close(&Zt, hzbuf);
      }

      dpd_buf4_mat_irrep_wrt(Z, hzbuf);
      dpd_buf4_mat_irrep_close(Z, hzbuf);

    }  /* end if(incore) */
    else { /* out-of-core for "normal" 424 contractions */
      /* Prepare the input buffer for the X factor and the target*/
#ifdef DPD_DEBUG	
      fprintf(stderr, "\t424 out-of-core: %d\n", hxbuf);
#endif
      dpd_buf4_mat_irrep_row_init(X, hxbuf);
      dpd_buf4_mat_irrep_row_init(Z, hzbuf);

      /* Loop over rows of the X factor and the target */
      for(pq=0; pq < Z->params->rowtot[hzbuf]; pq++) {

        dpd_buf4_mat_irrep_row_zero(X, hxbuf, pq);
        dpd_buf4_mat_irrep_row_rd(X, hxbuf, pq);

        dpd_buf4_mat_irrep_row_zero(Z, hzbuf, pq);

        if(fabs(beta) > 0.0)
          dpd_buf4_mat_irrep_row_rd(Z, hzbuf, pq);

        /*
           for(rs=0; rs < X->params->coltot[hxbuf]; rs++)
           fprintf(stderr, "\t%d X[%d] = %20.15f\n", pq, rs,
           X->matrix[hxbuf][0][rs]);


           for(rs=0; rs < Z->params->coltot[hzbuf]; rs++)
           fprintf(stderr, "\t%d Z[%d] = %20.15f\n", pq, rs,
           Z->matrix[hzbuf][0][rs]);
           fflush(stderr);
         */

        xcount = zcount = 0;

        for(Gr=0; Gr < nirreps; Gr++) {
          Gs = Gr^hxbuf;

          rowx = X->params->rpi[Gr];
          colx = X->params->spi[Gs];
          rowz = Z->params->rpi[Gr];
          colz = Z->params->spi[Gs];

          if(rowx && colx && colz) {
	    C_DGEMM('n',Ytrans?'y':'n',rowx,colz,colx,alpha,
		&(X->matrix[hxbuf][0][xcount]),colx,
		&(Y->matrix[Gs][0][0]),Ytrans?colx:colz,1.0,
		&(Z->matrix[hzbuf][0][zcount]),colz);
	  }
	  /*
            newmm2(&(X->matrix[hxbuf][0][xcount]),0,
                &(Y->matrix[Gs][0][0]),Ytrans,
                &(Z->matrix[hzbuf][0][zcount]),
                rowx,colx,colz,alpha,1.0);
		*/

          xcount += rowx * colx;
          zcount += rowz * colz;

        }

        dpd_buf4_mat_irrep_row_wrt(Z, hzbuf, pq);

        /*
           for(rs=0; rs < Z->params->coltot[hzbuf^GZ]; rs++)
           fprintf(stderr, "\t%d Zout[%d] = %20.15f\n", pq, rs,
           Z->matrix[hzbuf][0][rs]);
         */

      }

      dpd_buf4_mat_irrep_row_close(X, hxbuf);
      dpd_buf4_mat_irrep_row_close(Z, hzbuf);
    }
  }

  if((sum_X == 1) || (sum_X == 2)) dpd_trans4_close(&Xt);

  if(Ztrans) dpd_trans4_close(&Zt);

  dpd_file2_mat_close(Y);

  return 0;


}
