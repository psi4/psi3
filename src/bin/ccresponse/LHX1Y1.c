#include <stdio.h>
#include <string.h>
#include <libdpd/dpd.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#define EXTERN
#include "globals.h"

void build_ZIjAb(char *, char *, int, double, char *, char *, int, double);

double LHX1Y1(char *pert_x, char *cart_x, int irrep_x, double omega_x, 
	      char *pert_y, char *cart_y, int irrep_y, double omega_y)
{

  dpdfile2 F, X1, Y1, Zmi, Zae_1, Zae_2, Zfb, Znj, ZIA, L1, t1;
  dpdbuf4 Z1, Z2, I, tau, W1, W2, ZIjAb, L2, T2, W, Z;
  double polar;
  char lbl[32];
  int Gbm, Gfe, bm, b, m, Gb, Gm, Ge, Gf, B, M, fe, f, e, ef, nrows, ncols;
  double *X;

  /* The Lambda 1 contractions */
  dpd_file2_init(&ZIA, CC_TMP0, 0, 0, 1, "ZIA");
  sprintf(lbl, "X_%s_%1s_IA (%5.3f)", pert_x, cart_x, omega_x);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);
  sprintf(lbl, "X_%s_%1s_IA (%5.3f)", pert_y, cart_y, omega_y);
  dpd_file2_init(&Y1, CC_OEI, irrep_y, 0, 1, lbl);

  /* Contraction of FME, XIE, YMA */
  dpd_file2_init(&F, CC_OEI, 0, 0, 1, "FME");
  sprintf(lbl, "Z_%s_%1s_MI" , pert_x, cart_x);
  dpd_file2_init(&Zmi, CC_TMP0, irrep_x, 0, 0, lbl);
  dpd_contract222(&F, &X1, &Zmi, 0, 0, 1, 0);
  dpd_file2_close(&F);
  dpd_contract222(&Zmi, &Y1, &ZIA, 1, 1, -1, 0);

  /* Contraction of FME, XMA, YIE */
  dpd_file2_init(&F, CC_OEI, 0, 0, 1, "FME");
  dpd_file2_init(&Zmi, CC_TMP0, irrep_x, 0, 0, lbl);
  dpd_contract222(&F, &Y1, &Zmi, 0, 0, 1, 0);
  dpd_file2_close(&F);
  dpd_contract222(&Zmi, &X1, &ZIA, 1, 1, -1, 1);
  dpd_file2_close(&Zmi);

  /* Contraction of WAMEF with XIE, YMF and XMF, YIE */
  sprintf(lbl, "Z_%s_%1s_AE 1" , pert_y, cart_y);
  dpd_file2_init(&Zae_1, CC_TMP0, irrep_y, 1, 1, lbl);
  dpd_file2_scm(&Zae_1, 0);
  dpd_file2_mat_init(&Zae_1);
  sprintf(lbl, "Z_%s_%1s_AE 2" , pert_y, cart_y);
  dpd_file2_init(&Zae_2, CC_TMP0, irrep_y, 1, 1, lbl);
  dpd_file2_scm(&Zae_2, 0);
  dpd_file2_mat_init(&Zae_2);
  dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
  dpd_file2_mat_init(&Y1);
  dpd_file2_mat_rd(&Y1);
  dpd_file2_mat_init(&X1);
  dpd_file2_mat_rd(&X1);
  for(Gbm=0; Gbm < moinfo.nirreps; Gbm++) {
    Gfe = Gbm; /* W is totally symmetric */
    dpd_buf4_mat_irrep_row_init(&W, Gbm);
    X = init_array(W.params->coltot[Gfe]);
    for(bm=0; bm < W.params->rowtot[Gbm]; bm++) {
      dpd_buf4_mat_irrep_row_rd(&W, Gbm, bm);
      b = W.params->roworb[Gbm][bm][0];
      m = W.params->roworb[Gbm][bm][1];
      Gb = W.params->psym[b];
      Gm = Gbm ^ Gb;
      Ge = Gm ^ irrep_y;
      Gf = Gfe ^ Ge;
      B = b - moinfo.vir_off[Gb];
      M = m - moinfo.occ_off[Gm];
      zero_arr(X, W.params->coltot[Gfe]);
      for(fe=0; fe < W.params->coltot[Gfe]; fe++) {
	f = W.params->colorb[Gfe][fe][0];
	e = W.params->colorb[Gfe][fe][1];
	ef = W.params->colidx[e][f];
	X[fe] = 2.0 * W.matrix[Gbm][0][fe] - W.matrix[Gbm][0][ef];
      }
      nrows = moinfo.virtpi[Gf];
      ncols = moinfo.virtpi[Ge];
      if(nrows & ncols) {
	C_DGEMV('n',nrows,ncols,1,&X[W.col_offset[Gfe][Gf]],ncols,
		Y1.matrix[Gm][M],1,1,Zae_1.matrix[Gb][B],1);
	C_DGEMV('n',nrows,ncols,1,&X[W.col_offset[Gfe][Gf]],ncols,
		X1.matrix[Gm][M],1,1,Zae_2.matrix[Gb][B],1);
      }
    }
    free(X);
    dpd_buf4_mat_irrep_row_close(&W, Gbm);
  }
  dpd_file2_mat_close(&X1);
  dpd_file2_mat_close(&Y1);
  dpd_file2_mat_wrt(&Zae_1);
  dpd_file2_mat_close(&Zae_1);
  dpd_file2_mat_wrt(&Zae_2);
  dpd_file2_mat_close(&Zae_2);
  dpd_buf4_close(&W);
  dpd_contract222(&X1, &Zae_1, &ZIA, 0, 0, 1, 1);
  dpd_contract222(&Y1, &Zae_2, &ZIA, 0, 0, 1, 1);
  dpd_file2_close(&Zae_1);
  dpd_file2_close(&Zae_2);

  /* Contraction of WMNIE, XMA, YNE */
  sprintf(lbl, "Z_%s_%1s_MI" , pert_y, cart_y);
  dpd_file2_init(&Zmi, CC_TMP0, irrep_y, 0, 0, lbl);
  dpd_buf4_init(&W1, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe - 2WnMIe (Mn,eI)");
  dpd_dot13(&Y1, &W1, &Zmi, 0, 0, 1, 0);
  dpd_buf4_close(&W1);
  dpd_contract222(&Zmi, &X1, &ZIA, 1, 1, 1, 1);

  /* Contraction of WMNIE, XMA, YNE */
  dpd_buf4_init(&W1, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe - 2WnMIe (Mn,eI)");
  dpd_dot13(&X1, &W1, &Zmi, 0, 0, 1, 0);
  dpd_buf4_close(&W1);
  dpd_contract222(&Zmi, &Y1, &ZIA, 1, 1, 1, 1);
  dpd_file2_close(&Zmi);

  dpd_file2_close(&Y1);
  dpd_file2_close(&X1);

  /* Final contraction of ZIA intermediate with LIA */
  dpd_file2_init(&L1, CC_LAMPS, 0, 0, 1, "LIA 0 -1");
  polar = 2.0 * dpd_file2_dot(&ZIA, &L1);
  dpd_file2_close(&L1);
  dpd_file2_close(&ZIA);

  /*  fprintf(outfile, "L(1)HX1Y1 = %20.12f\n", polar); */

  /* The Lambda 2 contractions */
  sprintf(lbl, "X_%s_%1s_IA (%5.3f)", pert_x, cart_x, omega_x);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);
  sprintf(lbl, "X_%s_%1s_IA (%5.3f)", pert_y, cart_y, omega_y);
  dpd_file2_init(&Y1, CC_OEI, irrep_y, 0, 1, lbl);


  /* Contraction of Wmnij with Zmnab */
  dpd_buf4_init(&ZIjAb, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");
  dpd_buf4_scm(&ZIjAb, 0);
  build_ZIjAb(pert_x, cart_x, irrep_x, omega_x, pert_y, cart_y, irrep_y, omega_y);
  dpd_buf4_init(&Z1, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) anti");
  dpd_buf4_init(&W1, CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
  dpd_contract444(&W1, &Z1, &ZIjAb, 1, 1, 1, 0);
  dpd_buf4_close(&W1);
  dpd_buf4_close(&Z1);
  dpd_buf4_close(&ZIjAb);

  /* Contraction of Wabef with Zijef */
  dpd_buf4_init(&Z1, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) anti");

  dpd_buf4_init(&Z2, CC_TMP1, 0, 5, 0, 5, 0, 0, "Z(Ab,Ij)");
  dpd_buf4_init(&I, CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
  dpd_contract444(&I, &Z1, &Z2, 0, 0, 1, 0);
  dpd_buf4_close(&I);
  dpd_buf4_sort_axpy(&Z2, CC_TMP0, rspq, 0, 5, "Z(Ij,Ab) Final", 1);
  dpd_buf4_close(&Z2);

  dpd_buf4_init(&ZIjAb, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");
  dpd_buf4_init(&Z2, CC_TMP0, 0, 10, 0, 10, 0, 0, "Z(Mb,Ij)");
  dpd_buf4_init(&I, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_contract444(&I, &Z1, &Z2, 0, 0, -2, 0);
  dpd_buf4_close(&I);
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract244(&t1, &Z2, &ZIjAb, 0, 0, 1, 1, 1);
  dpd_file2_close(&t1);
  dpd_buf4_close(&Z2);
  dpd_buf4_init(&Z2, CC_TMP0, 0, 0, 0, 0, 0, 0, "Z(Ij,Mn)");
  dpd_buf4_init(&I, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_contract444(&Z1, &I, &Z2, 0, 0, 1, 0);
  dpd_buf4_close(&I);
  dpd_buf4_init(&tau, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
  dpd_contract444(&Z2, &tau, &ZIjAb, 0, 1, 1, 1);
  dpd_buf4_close(&tau);
  dpd_buf4_close(&Z2);
  dpd_buf4_close(&Z1);
  dpd_buf4_close(&ZIjAb);

  /* Contraction of Wmbej with Xie, Yma and Xma, Yie */
  dpd_buf4_init(&Z1, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) anti");
  dpd_buf4_sort(&Z1, CC_TMP0, psqr, 10, 10, "Z(Ib,jA) anti");
  dpd_buf4_close(&Z1);

  dpd_buf4_init(&Z2, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(IA,jb) I");
  dpd_buf4_init(&Z1, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(Ib,jA) anti");
  dpd_buf4_init(&W1, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
  dpd_contract444(&Z1, &W1, &Z2, 0, 1, -1, 0);
  dpd_buf4_close(&W1);
  dpd_buf4_close(&Z1);
  dpd_buf4_sort(&Z2, CC_TMP0, prqs, 0, 5, "Z(Ij,Ab) I");
  dpd_buf4_close(&Z2);

  dpd_buf4_init(&Z2, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(jA,Ib) II");
  dpd_buf4_init(&Z1, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(Ib,jA) anti");
  dpd_buf4_init(&W1, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
  dpd_contract444(&Z1, &W1, &Z2, 0, 1, 1, 0);
  dpd_buf4_close(&W1);
  dpd_buf4_close(&Z1);
  dpd_buf4_sort(&Z2, CC_TMP0, rpqs, 0, 5, "Z(Ij,Ab) II");
  dpd_buf4_close(&Z2);

  dpd_buf4_init(&Z1, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) I");
  dpd_buf4_init(&Z2, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) II");
  dpd_buf4_axpy(&Z2, &Z1, 1);
  dpd_buf4_close(&Z2);
  dpd_buf4_sort_axpy(&Z1, CC_TMP0, qpsr, 0, 5, "Z(Ij,Ab) Final", 1);
  dpd_buf4_init(&ZIjAb, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");
  dpd_buf4_axpy(&Z1, &ZIjAb, 1);
  dpd_buf4_close(&Z1);
  dpd_buf4_close(&ZIjAb);


  dpd_buf4_init(&ZIjAb, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");

  /* Contraction of Wmaneij with X and Y */
  dpd_file2_init(&Zfb, CC_TMP0, 0, 1, 1, "Z_fb");
  dpd_buf4_init(&Z1, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) anti");
  dpd_buf4_init(&I, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  dpd_contract442(&I, &Z1, &Zfb, 3, 3, -1, 0);
  dpd_buf4_close(&I);
  dpd_buf4_close(&Z1);

  dpd_buf4_init(&Z1, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) temp");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract424(&T2, &Zfb, &Z1, 3, 0, 0, 1, 0);
  dpd_buf4_close(&T2);
  dpd_file2_close(&Zfb);
  dpd_buf4_axpy(&Z1, &ZIjAb, 1);
  dpd_buf4_sort(&Z1, CC_TMP0, qpsr, 0, 5, "Z(jI,bA) temp");
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z1, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(jI,bA) temp");
  dpd_buf4_axpy(&Z1, &ZIjAb, 1);
  dpd_buf4_close(&Z1);

  /* Contraction of Wmabeif with X and Y */
  dpd_file2_init(&Znj, CC_TMP0, 0, 0, 0, "Z_nj");
  dpd_buf4_init(&Z1, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) anti");
  dpd_buf4_init(&I, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  dpd_contract442(&I, &Z1, &Znj, 1, 1, 1, 0);
  dpd_buf4_close(&I);
  dpd_buf4_close(&Z1);

  dpd_buf4_init(&Z1, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) temp");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract244(&Znj, &T2, &Z1, 0, 0, 0, -1, 0);
  dpd_buf4_close(&T2);
  dpd_file2_close(&Znj);
  dpd_buf4_axpy(&Z1, &ZIjAb, 1);
  dpd_buf4_sort(&Z1, CC_TMP0, qpsr, 0, 5, "Z(jI,bA) temp");
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z1, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(jI,bA) temp");
  dpd_buf4_axpy(&Z1, &ZIjAb, 1);
  dpd_buf4_close(&Z1);

  /* Close the X and Y matices */
  dpd_file2_close(&Y1);
  dpd_file2_close(&X1);

  /* Final contraction with LIJAB */
  dpd_buf4_init(&L2, CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
  polar += dpd_buf4_dot(&L2, &ZIjAb);
  dpd_buf4_close(&L2);
  dpd_buf4_close(&ZIjAb);

  return polar;
}
