/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libqt/qt.h>
#include "dpd.h"

extern "C" {

extern FILE *outfile;

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
 **
 ** applicability update by RAK, 2-04
 ** this function does not work for dots that require a 13 _and_ a 31 shift
 ** So currently we are using dots (0,0) (1,1) (2,2) (3,3) and (0,2);
 ** but (3,2) and (0,3) will not work unless C1 symmetry
 */

int dpd_contract442(dpdbuf4 *X, dpdbuf4 *Y, dpdfile2 *Z, int target_X,
    int target_Y, double alpha, double beta)
{
  int h,hxbuf,hybuf,nirreps,Gtar,GX,GY,GZ,Hx,Hy,Hz,hlinks;
  int rking=0;
  dpdtrans4 Xt, Yt;
  double ***Xmat, ***Ymat, ***Zmat;
  int Xtrans, Ytrans, *numlinks;
#ifdef DPD_DEBUG
  int *xrow, *xcol, *yrow, *ycol, *zrow, *zcol;
#endif

  nirreps = X->params->nirreps;
  GX = X->file.my_irrep;
  GY = Y->file.my_irrep;
  GZ = Z->my_irrep;

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

  /* loop over row buffer irreps of X */
  for(hxbuf=0; hxbuf < nirreps; hxbuf++) {
    if(target_X == 0) {
      dpd_buf4_mat_irrep_init(X, hxbuf);
      dpd_buf4_mat_irrep_rd(X, hxbuf);
      dpd_buf4_mat_irrep_shift13(X, hxbuf);
      Xmat = X->shift.matrix[hxbuf];
      Xtrans = 0;
      numlinks = X->shift.coltot[hxbuf];
#ifdef DPD_DEBUG
      xrow = X->shift.rowtot[hxbuf];
      xcol = X->shift.coltot[hxbuf];
#endif
    }
    else if(target_X == 1) {
      dpd_buf4_mat_irrep_init(X, hxbuf);
      dpd_buf4_mat_irrep_rd(X, hxbuf);
      dpd_trans4_mat_irrep_init(&Xt, hxbuf);
      dpd_trans4_mat_irrep_rd(&Xt, hxbuf);
      dpd_buf4_mat_irrep_close(X, hxbuf);
      dpd_trans4_mat_irrep_shift31(&Xt, hxbuf);
      rking = 1;
      Xmat = Xt.shift.matrix[hxbuf];
      Xtrans = 1;
      numlinks = Xt.shift.rowtot[hxbuf];
#ifdef DPD_DEBUG	  
      xrow = Xt.shift.coltot[hxbuf];
      xcol = Xt.shift.rowtot[hxbuf];
#endif
    }
    else if(target_X == 2) {
      dpd_buf4_mat_irrep_init(X, hxbuf);
      dpd_buf4_mat_irrep_rd(X, hxbuf);
      dpd_trans4_mat_irrep_init(&Xt, hxbuf);
      dpd_trans4_mat_irrep_rd(&Xt, hxbuf);
      dpd_buf4_mat_irrep_close(X, hxbuf);
      dpd_trans4_mat_irrep_shift13(&Xt, hxbuf);
      Xmat = Xt.shift.matrix[hxbuf];
      Xtrans = 0;
      numlinks = Xt.shift.coltot[hxbuf];
#ifdef DPD_DEBUG	  
      xrow = Xt.shift.rowtot[hxbuf];
      xcol = Xt.shift.coltot[hxbuf];
#endif	  
    }
    else if(target_X == 3) {
      dpd_buf4_mat_irrep_init(X, hxbuf);
      dpd_buf4_mat_irrep_rd(X, hxbuf);
      dpd_buf4_mat_irrep_shift31(X, hxbuf);
      rking = 1;
      Xmat = X->shift.matrix[hxbuf];
      Xtrans = 1;
      numlinks = X->shift.rowtot[hxbuf];
#ifdef DPD_DEBUG	  
      xrow = X->shift.coltot[hxbuf];
      xcol = X->shift.rowtot[hxbuf];
#endif	    
    }
    else {
      fprintf(stderr, "Junk X index %d in dpd_contract442\n", target_X);
      exit(PSI_RETURN_FAILURE);
    }

    /* read in appropriate block of Y buffer */
    if (target_X < 2) {
      if (target_Y < 2) hybuf = hxbuf^GZ; else hybuf = hxbuf^GX;
    }
    else {
      if (target_Y < 2) hybuf = hxbuf^GY; else hybuf = hxbuf;
    }

    if(target_Y == 0) {
      dpd_buf4_mat_irrep_init(Y, hybuf);
      dpd_buf4_mat_irrep_rd(Y, hybuf);
      dpd_buf4_mat_irrep_shift13(Y, hybuf);
      Ymat = Y->shift.matrix[hybuf];
      Ytrans = 1;
#ifdef DPD_DEBUG	  
      yrow = Y->shift.coltot[hybuf];
      ycol = Y->shift.rowtot[hybuf];
#endif		  
    }
    else if(target_Y == 1) {
      dpd_buf4_mat_irrep_init(Y, hybuf);
      dpd_buf4_mat_irrep_rd(Y, hybuf);
      dpd_trans4_mat_irrep_init(&Yt, hybuf);
      dpd_trans4_mat_irrep_rd(&Yt, hybuf);
      dpd_buf4_mat_irrep_close(Y, hybuf);
      dpd_trans4_mat_irrep_shift31(&Yt, hybuf);
      rking = 1;
      Ymat = Yt.shift.matrix[hybuf];
      Ytrans = 0;
#ifdef DPD_DEBUG	  
      yrow = Yt.shift.rowtot[hybuf];
      ycol = Yt.shift.coltot[hybuf];
#endif	  
    }
    else if(target_Y == 2) {
      dpd_buf4_mat_irrep_init(Y, hybuf);
      dpd_buf4_mat_irrep_rd(Y, hybuf);
      dpd_trans4_mat_irrep_init(&Yt, hybuf);
      dpd_trans4_mat_irrep_rd(&Yt, hybuf);
      dpd_buf4_mat_irrep_close(Y, hybuf);
      dpd_trans4_mat_irrep_shift13(&Yt, hybuf);
      Ymat = Yt.shift.matrix[hybuf];
      Ytrans = 1;
#ifdef DPD_DEBUG	  
      yrow = Yt.shift.coltot[hybuf];
      ycol = Yt.shift.rowtot[hybuf];
#endif	
    }
    else if(target_Y == 3) {
      dpd_buf4_mat_irrep_init(Y, hybuf);
      dpd_buf4_mat_irrep_rd(Y, hybuf);
      dpd_buf4_mat_irrep_shift31(Y, hybuf);
      rking = 1;
      Ymat = Y->shift.matrix[hybuf];
      Ytrans = 0;
#ifdef DPD_DEBUG	  
      yrow = Y->shift.rowtot[hybuf];
      ycol = Y->shift.coltot[hybuf];
#endif	  
    }
    else {
      fprintf(stderr, "Junk Y index %d in contract442\n", target_Y);
      exit(PSI_RETURN_FAILURE);
    }

    if(rking)
      for(Hx=0; Hx < nirreps; Hx++) {
#ifdef DPD_DEBUG
        if((xrow[Hx] != zrow[Hx]) || (ycol[Hx] != zcol[Hx]) ||
            (xcol[Hx] != yrow[Hx])) {
          fprintf(stderr, "** Alignment error in contract442 **\n");
          fprintf(stderr, "** Irrep %d; Subirrep %d **\n", h,Hx);
          dpd_error("dpd_contract442", stderr);
        }
#endif
        if      ((!Xtrans) && (!Ytrans)) { Hy = Hx^GX;    Hz = Hx;    }
        else if ((!Xtrans) && (Ytrans) ) { Hy = Hx^GX^GY; Hz = Hx;    }
        else if ( (Xtrans) && (!Ytrans)) { Hy = Hx;       Hz = Hx^GX; }
        else /* ( (Xtrans) && (Ytrans))*/{ Hy = Hx^GY;    Hz = Hx^GX; }

        /*
        if (!Xtrans && !Ytrans) {
	fprintf(outfile,"rows[%d]=%d links[%d]=%d cols[%d]=%d\n",
	   Hz, Z->params->rowtot[Hz], Hx, numlinks[Hx], Hz^GZ,
       Z->params->coltot[Hz^GZ]); 
    fflush(outfile);
    print_mat(Ymat[Hy], numlinks[Hx], Z->params->coltot[Hz^GZ], outfile);
    fflush(outfile);
        }
        */

        newmm_rking(Xmat[Hx], Xtrans, Ymat[Hy], Ytrans,
            Z->matrix[Hz], Z->params->rowtot[Hz],
            numlinks[Hx], Z->params->coltot[Hz^GZ],
            alpha, 1.0);
      }
    else
      for(Hx=0; Hx < nirreps; Hx++) {
#ifdef DPD_DEBUG
        if((xrow[Hx] != zrow[Hx]) || (ycol[Hx] != zcol[Hx]) ||
            (xcol[Hx] != yrow[Hx])) {
          fprintf(stderr, "** Alignment error in contract442 **\n");
          fprintf(stderr, "** Irrep %d; Subirrep %d **\n", h,Hx);
          dpd_error("dpd_contract442", stderr);
        }
#endif	      
        if      ((!Xtrans) && (!Ytrans)) { Hy = Hx^GX;    Hz = Hx;    }
        else if ((!Xtrans) && (Ytrans) ) { Hy = Hx^GX^GY; Hz = Hx;    }
        else if ( (Xtrans) && (!Ytrans)) { Hy = Hx;       Hz = Hx^GX; }
        else /* ( (Xtrans) && (Ytrans))*/{ Hy = Hx^GY;    Hz = Hx^GX; }
	/* fprintf(stdout,"rows %d links %d cols %d\n",
	   Z->params->rowtot[Hz], numlinks[Hx], Z->params->coltot[Hz]); */

	if(Z->params->rowtot[Hz] && 
	   Z->params->coltot[Hz^GZ] && 
	   numlinks[Hx]) {
	if(!Xtrans && !Ytrans) {
	  C_DGEMM('n','n',Z->params->rowtot[Hz],Z->params->coltot[Hz^GZ],
	      numlinks[Hx],alpha,&(Xmat[Hx][0][0]),numlinks[Hx],
	      &(Ymat[Hy][0][0]),Z->params->coltot[Hz^GZ],1.0,
	      &(Z->matrix[Hz][0][0]),Z->params->coltot[Hz^GZ]);
	}
	else if(Xtrans && !Ytrans) {
	  C_DGEMM('t','n',Z->params->rowtot[Hz],Z->params->coltot[Hz^GZ],
	      numlinks[Hx],alpha,&(Xmat[Hx][0][0]),Z->params->rowtot[Hz],
	      &(Ymat[Hy][0][0]),Z->params->coltot[Hz^GZ],1.0,
	      &(Z->matrix[Hz][0][0]),Z->params->coltot[Hz^GZ]);
	}
	else if(!Xtrans && Ytrans) {
	  C_DGEMM('n','t',Z->params->rowtot[Hz],Z->params->coltot[Hz^GZ],
	      numlinks[Hx],alpha,&(Xmat[Hx][0][0]),numlinks[Hx],
	      &(Ymat[Hy][0][0]),numlinks[Hx],1.0,
	      &(Z->matrix[Hz][0][0]),Z->params->coltot[Hz^GZ]);
	}
	else {
	  C_DGEMM('t','t',Z->params->rowtot[Hz],Z->params->coltot[Hz^GZ],
	      numlinks[Hx],alpha,&(Xmat[Hx][0][0]),Z->params->rowtot[Hz],
	      &(Ymat[Hy][0][0]),numlinks[Hx],1.0,
	      &(Z->matrix[Hz][0][0]),Z->params->coltot[Hz^GZ]);
	}
	}
	/*
        newmm(Xmat[Hx], Xtrans, Ymat[Hy], Ytrans,
            Z->matrix[Hz], Z->params->rowtot[Hz],
            numlinks[Hx], Z->params->coltot[Hz^GZ],
            alpha, 1.0);
	    */
      }

    if(target_X == 0) dpd_buf4_mat_irrep_close(X, hxbuf);
    else if(target_X == 1) dpd_trans4_mat_irrep_close(&Xt, hxbuf);
    else if(target_X == 2) dpd_trans4_mat_irrep_close(&Xt, hxbuf);
    else if(target_X == 3) dpd_buf4_mat_irrep_close(X, hxbuf);

    if(target_Y == 0) dpd_buf4_mat_irrep_close(Y, hybuf);
    else if(target_Y == 1) dpd_trans4_mat_irrep_close(&Yt, hybuf);
    else if(target_Y == 2) dpd_trans4_mat_irrep_close(&Yt, hybuf);
    else if(target_Y == 3) dpd_buf4_mat_irrep_close(Y, hybuf);
  }

  if((target_X == 1) || (target_X == 2)) dpd_trans4_close(&Xt);
  if((target_Y == 1) || (target_Y == 2)) dpd_trans4_close(&Yt);


  dpd_file2_mat_wrt(Z);
  dpd_file2_mat_close(Z);

  return 0;
}



} /* extern "C" */
