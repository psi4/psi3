#include <stdio.h>
#include <math.h>
#include <libqt/qt.h>
#include "dpd.h"
#define EXTERN
#include "dpd.gbl"
extern FILE *outfile;

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
 **   int Ztrans: A boolean which indicates whether the indices in Z
 **                must be transposed in order to match the external
 **                index ordering in the product X * Y.
 **   double alpha: A prefactor for the product alpha * X * Y.
 **   double beta: A prefactor for the target beta * Z.
 */

int dpd_contract244(dpdfile2 *X, dpdbuf4 *Y, dpdbuf4 *Z, int sum_X, int sum_Y,
    int Ztrans, double alpha, double beta)
{
  int h, h0, hx, hybuf, hzbuf, hy, hz, nirreps, GX, GY, GZ, bra_y;
  int rking=0, *yrow, *ycol, symlink;
  int Xtrans, Ytrans = 0;
  int *numlinks, *numrows, *numcols;
  dpdtrans4 Yt, Zt;
  double ***Ymat, ***Zmat;
#ifdef DPD_DEBUG
  int *xrow, *xcol, *zrow, *zcol;
#endif  

  nirreps = X->params->nirreps;
  GX = X->my_irrep;
  GY = Y->file.my_irrep;
  GZ = Z->file.my_irrep;

  dpd_file2_mat_init(X);
  dpd_file2_mat_rd(X);

  if(sum_X == 0) { Xtrans = 1; numlinks = X->params->rowtot; symlink = 0; }
  else if(sum_X == 1) { Xtrans = 0; numlinks = X->params->coltot; symlink = GX; }
  else { fprintf(stderr, "Junk X index %d\n", sum_X); exit(sum_X); }

  if((sum_Y == 1) || (sum_Y == 2)) dpd_trans4_init(&Yt, Y);

  if(Ztrans) dpd_trans4_init(&Zt, Z);

  /*  if(fabs(beta) > 0.0) dpd_buf4_scm(Z, beta); */
  dpd_buf4_scm(Z, beta);

#ifdef DPD_DEBUG  
  if(Xtrans) { xrow = X->params->coltot; xcol = X->params->rowtot; }
  else { xrow = X->params->rowtot; xcol = X->params->coltot; }
#endif  

  for(hzbuf=0; hzbuf < nirreps; hzbuf++) {

    dpd_buf4_mat_irrep_init(Z, hzbuf);
    if(fabs(beta) > 0.0) dpd_buf4_mat_irrep_rd(Z, hzbuf);
    if(Ztrans) {
      dpd_trans4_mat_irrep_init(&Zt, hzbuf);
      dpd_trans4_mat_irrep_rd(&Zt, hzbuf);
      dpd_buf4_mat_irrep_close(Z, hzbuf);
      dpd_trans4_mat_irrep_shift13(&Zt, hzbuf);
      numrows = Zt.shift.rowtot[hzbuf];
      numcols = Zt.shift.coltot[hzbuf];
      Zmat = Zt.shift.matrix[hzbuf];
#ifdef DPD_DEBUG
      zrow = Zt.shift.rowtot[hzbuf];
      zcol = Zt.shift.coltot[hzbuf];
#endif
    }
    else {
      dpd_buf4_mat_irrep_shift13(Z, hzbuf);
      numrows = Z->shift.rowtot[hzbuf];
      numcols = Z->shift.coltot[hzbuf];
      Zmat = Z->shift.matrix[hzbuf];
#ifdef DPD_DEBUG
      zrow = Z->shift.rowtot[hzbuf];
      zcol = Z->shift.coltot[hzbuf];
#endif	  
    }

    if (sum_Y < 2) {
      if (Ztrans) hybuf = hzbuf^GY; else hybuf = hzbuf^GZ^GY;
    }
    else {
      if (Ztrans) hybuf = hzbuf; else hybuf = hzbuf^GZ;
    }
       
    if(sum_Y == 0) {
      dpd_buf4_mat_irrep_init(Y, hybuf);
      dpd_buf4_mat_irrep_rd(Y, hybuf);
      dpd_buf4_mat_irrep_shift13(Y, hybuf);
      Ymat = Y->shift.matrix[hybuf];
      Ytrans = 0;
#ifdef DPD_DEBUG
      yrow = Y->shift.rowtot[hybuf];
      ycol = Y->shift.coltot[hybuf];
#endif
    }
    else if(sum_Y == 1) {
      dpd_buf4_mat_irrep_init(Y, hybuf);
      dpd_buf4_mat_irrep_rd(Y, hybuf);
      dpd_trans4_mat_irrep_init(&Yt, hybuf);
      dpd_trans4_mat_irrep_rd(&Yt, hybuf);
      dpd_buf4_mat_irrep_close(Y, hybuf);
      dpd_trans4_mat_irrep_shift31(&Yt, hybuf);
      rking = 1;
      Ytrans = 1;
      Ymat = Yt.shift.matrix[hybuf];
#ifdef DPD_DEBUG
      yrow = Yt.shift.coltot[hybuf];
      ycol = Yt.shift.rowtot[hybuf];
#endif  
    }
    else if(sum_Y == 2) {
      dpd_buf4_mat_irrep_init(Y, hybuf);
      dpd_buf4_mat_irrep_rd(Y, hybuf);
      dpd_trans4_mat_irrep_init(&Yt, hybuf);
      dpd_trans4_mat_irrep_rd(&Yt, hybuf);
      dpd_buf4_mat_irrep_close(Y, hybuf);
      dpd_trans4_mat_irrep_shift13(&Yt, hybuf);
      Ymat = Yt.shift.matrix[hybuf];
      Ytrans = 0;
#ifdef DPD_DEBUG
      yrow = Yt.shift.rowtot[hybuf];
      ycol = Yt.shift.coltot[hybuf];
#endif 	  
    }
    else if(sum_Y == 3) {
      dpd_buf4_mat_irrep_init(Y, hybuf);
      dpd_buf4_mat_irrep_rd(Y, hybuf);
      dpd_buf4_mat_irrep_shift31(Y, hybuf);
      rking = 1;
      Ytrans = 1;
      Ymat = Y->shift.matrix[hybuf];
#ifdef DPD_DEBUG
      yrow = Y->shift.coltot[hybuf];
      ycol = Y->shift.rowtot[hybuf];
#endif		  
    }

    if(rking) 
      for(hz=0; hz < nirreps; hz++) {
        if      (!Xtrans && !Ytrans) {hx=hz;       hy = hz^GX; }
        else if (!Xtrans &&  Ytrans) {hx=hz;       hy = hz^GX^GY; }
        else if ( Xtrans && !Ytrans) {hx=hz^GX;    hy = hz^GX; }
        else if ( Xtrans &&  Ytrans) {hx=hz^GX;    hy = hz^GX^GY; }
#ifdef DPD_DEBUG
        if((xrow[hz] != zrow[hz]) || (ycol[hz] != zcol[hz]) ||
            (xcol[hz] != yrow[hz])) {
          fprintf(stderr, "** Alignment error in contract244 **\n");
          fprintf(stderr, "** Irrep: %d; Subirrep: %d **\n",hzbuf,hz);
          dpd_error("dpd_contract244", stderr);
        }
#endif
        newmm_rking(X->matrix[hx],Xtrans, Ymat[hy], Ytrans,
            Zmat[hz], numrows[hz], numlinks[hx^symlink],
            numcols[hz], alpha, 1.0);
      }
    else
      for(hz=0; hz < nirreps; hz++) {
        if      (!Xtrans && !Ytrans ) {hx=hz;       hy = hz^GX; }
        else if (!Xtrans &&  Ytrans ) {hx=hz;       hy = hz^GX^GY; }
        else if ( Xtrans && !Ytrans ) {hx=hz^GX;    hy = hz^GX; }
        else if ( Xtrans &&  Ytrans ) {hx=hz^GX;    hy = hz^GX^GY; }

#ifdef DPD_DEBUG
        if((xrow[hz] != zrow[hz]) || (ycol[hz] != zcol[hz]) ||
            (xcol[hz] != yrow[hz])) {
          fprintf(stderr, "** Alignment error in contract244 **\n");
          fprintf(stderr, "** Irrep: %d; Subirrep: %d **\n",hzbuf,hz);
          dpd_error("dpd_contract244", stderr);
        }
#endif	      
	/* fprintf(outfile,"hz %d, hx %d, hy %d, numrows %d, numlinks %d, numcols %d\n",
	   hz, hx, hy, numrows[hz],numlinks[hx],numcols[hz]); */

        newmm(X->matrix[hx], Xtrans, Ymat[hy], Ytrans,
            Zmat[hz], numrows[hz], numlinks[hx^symlink],
            numcols[hz], alpha, 1.0);
      }

    if(sum_Y == 0) dpd_buf4_mat_irrep_close(Y, hybuf);
    else if(sum_Y == 1) dpd_trans4_mat_irrep_close(&Yt, hybuf);
    else if(sum_Y == 2) dpd_trans4_mat_irrep_close(&Yt, hybuf);
    else if(sum_Y == 3) dpd_buf4_mat_irrep_close(Y, hybuf);

    if(Ztrans) {
      dpd_buf4_mat_irrep_init(Z, hzbuf);
      dpd_trans4_mat_irrep_wrt(&Zt, hzbuf);
      dpd_trans4_mat_irrep_close(&Zt, hzbuf);
    }

    dpd_buf4_mat_irrep_wrt(Z, hzbuf);
    dpd_buf4_mat_irrep_close(Z, hzbuf);

  }

  if((sum_Y == 1) || (sum_Y == 2)) dpd_trans4_close(&Yt);

  if(Ztrans) dpd_trans4_close(&Zt);

  dpd_file2_mat_close(X);

  return 0;
}
