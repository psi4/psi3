#include <stdio.h>
#include <string.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void build_ZIjAb(char *, int, double);

double LHX1Y1(char *cart, int irrep, double omega) {

  dpdfile2 F, X1, Y1, Zmi, Zae, Zfb, Znj, ZIA, L1, t1;
  dpdbuf4 Z1, Z2, I, tau, W1, W2, ZIjAb, L2, T2, W, Z;
  double polar;
  char lbl[32];

  /* The Lambda 1 contractions */
  dpd_file2_init(&ZIA, CC_TMP0, 0, 0, 1, "ZIA");
  sprintf(lbl, "X_%1s_IA (%5.3f)", cart, omega);
  dpd_file2_init(&X1, CC_OEI, irrep, 0, 1, lbl);
  sprintf(lbl, "X_%1s_IA (-%5.3f)", cart, omega);
  dpd_file2_init(&Y1, CC_OEI, irrep, 0, 1, lbl);

  /* Contraction of FME, XIE, YMA */
  dpd_file2_init(&F, CC_OEI, 0, 0, 1, "FME");
  sprintf(lbl, "Z_%1s_MI" , cart);
  dpd_file2_init(&Zmi, CC_TMP0, irrep, 0, 0, lbl);
  dpd_contract222(&F, &X1, &Zmi, 0, 0, 1, 0);
  dpd_file2_close(&F);
  dpd_contract222(&Zmi, &Y1, &ZIA, 1, 1, -1, 0);

  /* Contraction of FME, XMA, YIE */
  dpd_file2_init(&F, CC_OEI, 0, 0, 1, "FME");
  dpd_file2_init(&Zmi, CC_TMP0, irrep, 0, 0, lbl);
  dpd_contract222(&F, &Y1, &Zmi, 0, 0, 1, 0);
  dpd_file2_close(&F);
  dpd_contract222(&Zmi, &X1, &ZIA, 1, 1, -1, 1);
  dpd_file2_close(&Zmi);

  /* Contraction of WAMEF, XIE, YMF */
  sprintf(lbl, "Z_%1s_AE" , cart);
  dpd_file2_init(&Zae, CC_TMP0, irrep, 1, 1, lbl);
  dpd_buf4_init(&W1, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)");
  dpd_dot24(&Y1, &W1, &Zae, 0, 0, 1, 0);
  dpd_buf4_close(&W1);
  dpd_contract222(&X1, &Zae, &ZIA, 0, 0, 1, 1);

  /* Contraction of WAMEF, XMF, YIE */
  dpd_buf4_init(&W1, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)");
  dpd_dot24(&X1, &W1, &Zae, 0, 0, 1, 0);
  dpd_buf4_close(&W1);
  dpd_contract222(&Y1, &Zae, &ZIA, 0, 0, 1, 1);
  dpd_file2_close(&Zae);

  /* Contraction of WAMEF, XMA, YNE */
  sprintf(lbl, "Z_%1s_MI" , cart);
  dpd_file2_init(&Zmi, CC_TMP0, irrep, 0, 0, lbl);
  dpd_buf4_init(&W1, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe - 2WnMIe");
  dpd_dot13(&Y1, &W1, &Zmi, 0, 0, 1, 0);
  dpd_buf4_close(&W1);
  dpd_contract222(&Zmi, &X1, &ZIA, 1, 1, 1, 1);

  /* Contraction of WAMEF, XMA, YNE */
  dpd_buf4_init(&W1, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe - 2WnMIe");
  dpd_dot13(&X1, &W1, &Zmi, 0, 0, 1, 0);
  dpd_buf4_close(&W1);
  dpd_contract222(&Zmi, &Y1, &ZIA, 1, 1, 1, 1);
  dpd_file2_close(&Zmi);

  dpd_file2_close(&Y1);
  dpd_file2_close(&X1);

  /* Final contraction of ZIA intermediate with LIA */
  dpd_file2_init(&L1, CC_OEI, 0, 0, 1, "LIA");
  polar = 2.0 * dpd_file2_dot(&ZIA, &L1);
  dpd_file2_close(&L1);
  dpd_file2_close(&ZIA);

  /*  fprintf(outfile, "L(1)HX1Y1 = %20.12f\n", polar); */

  /* The Lambda 2 contractions */
  sprintf(lbl, "X_%1s_IA (%5.3f)", cart, omega);
  dpd_file2_init(&X1, CC_OEI, irrep, 0, 1, lbl);
  sprintf(lbl, "X_%1s_IA (-%5.3f)", cart, omega);
  dpd_file2_init(&Y1, CC_OEI, irrep, 0, 1, lbl);


  /* Contraction of Wmnij with Zmnab */
  dpd_buf4_init(&ZIjAb, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");
  dpd_buf4_scm(&ZIjAb, 0);
  build_ZIjAb(cart, irrep, omega);
  dpd_buf4_init(&Z1, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) anti");
  dpd_buf4_init(&W1, CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
  dpd_contract444(&W1, &Z1, &ZIjAb, 1, 1, 1, 0);
  dpd_buf4_close(&W1);
  dpd_buf4_close(&Z1);
  dpd_buf4_close(&ZIjAb);

  /* Contraction of Wabef with Zijef */
  dpd_buf4_init(&ZIjAb, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");
  dpd_buf4_init(&Z1, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) anti");
  dpd_buf4_init(&I, CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
  dpd_contract444(&Z1, &I, &ZIjAb, 0, 0, 1, 1);
  dpd_buf4_close(&I);
  dpd_buf4_init(&Z2, CC_TMP0, 0, 0, 11, 0, 11, 0, "Z(ij,am)");
  dpd_buf4_init(&I, CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
  dpd_contract444(&Z1, &I, &Z2, 0, 0, -2, 0);
  dpd_buf4_close(&I);
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&Z2, &t1, &ZIjAb, 3, 0, 0, 1, 1);
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
