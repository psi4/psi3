#include <stdio.h>
#include <string.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

double LHX1Y2(char *cart_x, int irrep_x, double omega_x, char *cart_y, int irrep_y, double omega_y)
{
  dpdfile2 z, z1, X1, l1, F;
  dpdbuf4 Z, Z1, Z2, I, Y2, L2, W;
  char lbl[32];
  double polar;

  sprintf(lbl, "Z_%1s_MI", cart_y);
  dpd_file2_init(&z1, CC_TMP0, irrep_y, 0, 0, lbl);
  dpd_buf4_init(&I, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  sprintf(lbl, "X_%1s_IjAb (%5.3f)", cart_y, omega_y);
  dpd_buf4_init(&Y2, CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);
  dpd_contract442(&I, &Y2, &z1, 0, 0, 1, 0);
  dpd_buf4_close(&Y2);
  dpd_buf4_close(&I);

  dpd_file2_init(&z, CC_TMP0, 0, 0, 1, "Z(I,A) Final");
  sprintf(lbl, "X_%1s_IA (%5.3f)", cart_x, omega_x);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);
  dpd_contract222(&z1, &X1, &z, 1, 1, -1, 0);
  dpd_file2_close(&X1);
  dpd_file2_close(&z1);
  dpd_file2_close(&z);

  sprintf(lbl, "Z_%1s_AE", cart_y);
  dpd_file2_init(&z1, CC_TMP0, irrep_y, 1, 1, lbl);
  dpd_buf4_init(&I, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  sprintf(lbl, "X_%1s_IjAb (%5.3f)", cart_y, omega_y);
  dpd_buf4_init(&Y2, CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);
  dpd_contract442(&Y2, &I, &z1, 3, 3, -1, 0);
  dpd_buf4_close(&Y2);
  dpd_buf4_close(&I);

  dpd_file2_init(&z, CC_TMP0, 0, 0, 1, "Z(I,A) Final");
  sprintf(lbl, "X_%1s_IA (%5.3f)", cart_x, omega_x);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);
  dpd_contract222(&X1, &z1, &z, 0, 0, 1, 1);
  dpd_file2_close(&X1);
  dpd_file2_close(&z1);
  dpd_file2_close(&z);


  sprintf(lbl, "Z_%1s_ME", cart_x);
  dpd_file2_init(&z1, CC_TMP0, irrep_x, 0, 1, lbl);
  sprintf(lbl, "X_%1s_IA (%5.3f)", cart_x, omega_x);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);
  dpd_buf4_init(&I, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  dpd_dot24(&X1, &I, &z1, 0, 0, 1, 0);
  dpd_buf4_close(&I);
  dpd_file2_close(&X1);

  dpd_file2_init(&z, CC_TMP0, 0, 0, 1, "Z(I,A) Final");
  sprintf(lbl, "X_%1s_(2IjAb-IjbA) (%5.3f)", cart_y, omega_y);
  dpd_buf4_init(&Y2, CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);
  dpd_dot24(&z1, &Y2, &z, 0, 0, 1, 1);
  dpd_buf4_close(&Y2);
  dpd_file2_close(&z1);
  dpd_file2_close(&z);

  dpd_file2_init(&z, CC_TMP0, 0, 0, 1, "Z(I,A) Final");
  dpd_file2_init(&l1, CC_OEI, 0, 0, 1, "LIA");
  polar = 2.0 * dpd_file2_dot(&z, &l1);
  dpd_file2_close(&l1);
  dpd_file2_close(&z);

  /*  fprintf(outfile, "L(1)HX1Y2 = %20.12f\n", polar); */


  dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");
  dpd_buf4_scm(&Z, 0);

  sprintf(lbl, "Z_%1s_MI", cart_x);
  dpd_file2_init(&z1, CC_TMP0, irrep_x, 0, 0, lbl);
  sprintf(lbl, "X_%1s_IA (%5.3f)", cart_x, omega_x);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);
  dpd_file2_init(&F, CC_OEI, 0, 0, 1, "FME");
  dpd_contract222(&F, &X1, &z1, 0, 0, 1, 0);
  dpd_file2_close(&F);
  dpd_file2_close(&X1);

  dpd_buf4_init(&Z1, CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
  sprintf(lbl, "X_%1s_IjAb (%5.3f)", cart_y, omega_y);
  dpd_buf4_init(&Y2, CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);
  dpd_contract244(&z1, &Y2, &Z1, 0, 0, 0, 1, 0);
  dpd_buf4_close(&Y2);
  dpd_file2_close(&z1);
  dpd_buf4_axpy(&Z1, &Z, -1);
  dpd_buf4_sort(&Z1, CC_TMP1, qpsr, 0, 5, "Z(jI,bA)");
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z1, CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(jI,bA)");
  dpd_buf4_axpy(&Z1, &Z, -1);
  dpd_buf4_close(&Z1);

  dpd_buf4_close(&Z);


  dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");

  sprintf(lbl, "Z_%1s_AE", cart_x);
  dpd_file2_init(&z1, CC_TMP0, irrep_x, 1, 1, lbl);
  sprintf(lbl, "X_%1s_IA (%5.3f)", cart_x, omega_x);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);
  dpd_file2_init(&F, CC_OEI, 0, 0, 1, "FME");
  dpd_contract222(&X1, &F, &z1, 1, 1, -1, 0);
  dpd_file2_close(&F);
  dpd_file2_close(&X1);

  dpd_buf4_init(&Z1, CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
  sprintf(lbl, "X_%1s_IjAb (%5.3f)", cart_y, omega_y);
  dpd_buf4_init(&Y2, CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);
  dpd_contract424(&Y2, &z1, &Z1, 3, 1, 0, 1, 0);
  dpd_buf4_close(&Y2);
  dpd_file2_close(&z1);
  dpd_buf4_axpy(&Z1, &Z, 1);
  dpd_buf4_sort(&Z1, CC_TMP1, qpsr, 0, 5, "Z(jI,bA)");
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z1, CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(jI,bA)");
  dpd_buf4_axpy(&Z1, &Z, 1);
  dpd_buf4_close(&Z1);

  dpd_buf4_close(&Z);

  sprintf(lbl, "X_%1s_IA (%5.3f)", cart_x, omega_x);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);

  sprintf(lbl, "Z_%1s_MbEj (Mb,jE)", cart_x);
  dpd_buf4_init(&Z1, CC_TMP0, irrep_x, 10, 10, 10, 10, 0, lbl);
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 5, 10, 5, 0, "WAmEf");
  dpd_contract244(&X1, &W, &Z1, 1, 2, 1, 1, 0);
  dpd_buf4_close(&W);
  sprintf(lbl, "Z_%1s_MbEj (ME,jb)", cart_x);
  dpd_buf4_sort(&Z1, CC_TMP0, psrq, 10, 10, lbl);
  dpd_buf4_close(&Z1);

  sprintf(lbl, "Z_%1s_WMbEj (bM,Ej)", cart_x);
  dpd_buf4_init(&Z1, CC_TMP0, irrep_x, 11, 11, 11, 11, 0, lbl);
  dpd_buf4_init(&W, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe");
  dpd_contract244(&X1, &W, &Z1, 0, 0, 0, -1, 0);
  dpd_buf4_close(&W);
  sprintf(lbl, "Z_%1s_MbEj (ME,jb)", cart_x);
  dpd_buf4_sort_axpy(&Z1, CC_TMP0, qrsp, 10, 10, lbl, 1);
  dpd_buf4_close(&Z1);

  sprintf(lbl, "Z_%1s_MbeJ (Mb,eJ)", cart_x);
  dpd_buf4_init(&Z1, CC_TMP0, irrep_x, 10, 11, 10, 11, 0, lbl);
  dpd_buf4_init(&W, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe");
  dpd_contract424(&W, &X1, &Z1, 1, 0, 1, 1, 0);
  dpd_buf4_close(&W);
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 5, 10, 5, 0, "WAmEf");
  dpd_contract424(&W, &X1, &Z1, 3, 1, 0, -1, 1);
  dpd_buf4_close(&W);
  sprintf(lbl, "Z_%1s_MbeJ (Me,Jb)", cart_x);
  dpd_buf4_sort(&Z1, CC_TMP0, prsq, 10, 10, lbl);
  dpd_buf4_close(&Z1);

  dpd_file2_close(&X1);

  sprintf(lbl, "Z_%1s_MbEj (ME,jb)", cart_x);
  dpd_buf4_init(&Z1, CC_TMP0, irrep_x, 10, 10, 10, 10, 0, lbl);
  sprintf(lbl, "Z_%1s_(2MbEj+MbeJ) (ME,JB)", cart_x);
  dpd_buf4_scmcopy(&Z1, CC_TMP0, lbl, 2);
  dpd_buf4_close(&Z1);
  sprintf(lbl, "Z_%1s_MbeJ (Me,Jb)", cart_x);
  dpd_buf4_init(&Z1, CC_TMP0, irrep_x, 10, 10, 10, 10, 0, lbl);
  sprintf(lbl, "Z_%1s_(2MbEj+MbeJ) (ME,JB)", cart_x);
  dpd_buf4_init(&Z2, CC_TMP0, irrep_x, 10, 10, 10, 10, 0, lbl);
  dpd_buf4_axpy(&Z1, &Z2, 1);
  dpd_buf4_close(&Z2);
  dpd_buf4_close(&Z1);

  dpd_buf4_init(&Z1, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(IA,jb) I");
  sprintf(lbl, "X_%1s_(2IAjb-IbjA) (%5.3f)", cart_y, omega_y);
  dpd_buf4_init(&Y2, CC_LR, irrep_y, 10, 10, 10, 10, 0, lbl);
  sprintf(lbl, "Z_%1s_(2MbEj+MbeJ) (ME,JB)", cart_x);
  dpd_buf4_init(&Z, CC_TMP0, irrep_x, 10, 10, 10, 10, 0, lbl);
  dpd_contract444(&Y2, &Z, &Z1, 0, 1, 0.5, 0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&Y2);
  dpd_buf4_close(&Z1);

  dpd_buf4_init(&Z1, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(IA,jb) Ib");
  sprintf(lbl, "Z_%1s_MbeJ (Me,Jb)", cart_x);
  dpd_buf4_init(&Z, CC_TMP0, irrep_x, 10, 10, 10, 10, 0, lbl);
  sprintf(lbl, "X_%1s_IbjA (%5.3f)", cart_y, omega_y);
  dpd_buf4_init(&Y2, CC_LR, irrep_y, 10, 10, 10, 10, 0, lbl);
  dpd_contract444(&Y2, &Z, &Z1, 0, 1, 1, 0);
  dpd_buf4_close(&Y2);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&Z1);

  dpd_buf4_init(&Z1, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(IA,jb) I");
  dpd_buf4_init(&Z2, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(IA,jb) Ib");
  dpd_buf4_axpy(&Z2, &Z1, 0.5);
  dpd_buf4_sort(&Z2, CC_TMP0, psrq, 10, 10, "Z(IA,jb) III");
  dpd_buf4_close(&Z2);
  dpd_buf4_close(&Z1);

  dpd_buf4_init(&Z1, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(IA,jb) I");
  dpd_buf4_init(&Z2, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(IA,jb) III");
  dpd_buf4_axpy(&Z2, &Z1, 1);
  dpd_buf4_close(&Z2);
  dpd_buf4_sort(&Z1, CC_TMP0, prqs, 0, 5, "Z(Ij,Ab) I+III");
  dpd_buf4_close(&Z1);

  dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) I+III");
  dpd_buf4_sort(&Z, CC_TMP0, qpsr, 0, 5, "Z(Ij,Ab) II+IV");
  dpd_buf4_close(&Z);

  dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");
  dpd_buf4_init(&Z1, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) I+III");
  dpd_buf4_axpy(&Z1, &Z, 1);
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z1, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) II+IV");
  dpd_buf4_axpy(&Z1, &Z, 1);
  dpd_buf4_close(&Z1);
  dpd_buf4_close(&Z);


  sprintf(lbl, "X_%1s_IA (%5.3f)", cart_x, omega_x);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);
  sprintf(lbl, "X_%1s_IjAb (%5.3f)", cart_y, omega_y);
  dpd_buf4_init(&Y2, CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);

  sprintf(lbl, "Z_%1s_mj" , cart_x);
  dpd_file2_init(&z, CC_TMP0, irrep_x, 0, 0, lbl);
  dpd_buf4_init(&W, CC_HBAR, 0, 0, 11, 0, 11, 0, "2WMnIe - WnMIe");
  dpd_dot23(&X1, &W, &z, 0, 0, 1, 0);
  dpd_buf4_close(&W);
  dpd_buf4_init(&Z1, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
  dpd_contract424(&Y2, &z, &Z1, 1, 0, 1, -1, 0);
  dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");
  dpd_buf4_axpy(&Z1, &Z, 1);
  dpd_buf4_close(&Z);
  dpd_buf4_sort_axpy(&Z1, CC_TMP0, qpsr, 0, 5, "Z(Ij,Ab) Final", 1);
  dpd_buf4_close(&Z1);
  dpd_file2_close(&z);



  sprintf(lbl, "Z_%1s_ae" , cart_x);
  dpd_file2_init(&z, CC_TMP0, irrep_x, 1, 1, lbl);

  dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)");
  dpd_dot24(&X1, &W, &z, 0, 0, 1, 0);
  dpd_buf4_close(&W);
  dpd_buf4_init(&Z1, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) temp");
  dpd_contract424(&Y2, &z, &Z1, 3, 1, 0, 1, 0);
  dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");
  dpd_buf4_axpy(&Z1, &Z, 1);
  dpd_buf4_close(&Z);
  dpd_buf4_sort_axpy(&Z1, CC_TMP0, qpsr, 0, 5, "Z(Ij,Ab) Final", 1);
  dpd_buf4_close(&Z1);
  dpd_file2_close(&z);

  dpd_file2_close(&X1);
  dpd_buf4_close(&Y2);


  sprintf(lbl, "X_%1s_IA (%5.3f)", cart_x, omega_x);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);
  sprintf(lbl, "Z_%1s_MnjI", cart_x);
  dpd_buf4_init(&Z1, CC_TMP0, irrep_x, 0, 0, 0, 0, 0, lbl);
  dpd_buf4_init(&W, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe");
  dpd_contract244(&X1, &W, &Z1, 1, 2, 1, 1, 0);
  dpd_buf4_close(&W);
  sprintf(lbl, "Z_%1s_MnIj", cart_x);
  dpd_buf4_sort(&Z1, CC_TMP0, pqsr, 0, 0, lbl);
  dpd_buf4_sort_axpy(&Z1, CC_TMP0, qprs, 0, 0, lbl, 1);
  dpd_buf4_close(&Z1);
  dpd_file2_close(&X1);

  dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");
  sprintf(lbl, "Z_%1s_MnIj", cart_x);
  dpd_buf4_init(&Z1, CC_TMP0, irrep_x, 0, 0, 0, 0, 0, lbl);
  sprintf(lbl, "X_%1s_IjAb (%5.3f)", cart_y, omega_y);
  dpd_buf4_init(&Y2, CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);
  dpd_contract444(&Z1, &Y2, &Z, 1, 1, 1, 1);
  dpd_buf4_close(&Y2);
  dpd_buf4_close(&Z1);
  dpd_buf4_close(&Z);

  sprintf(lbl, "Z_%1s_AmIj (mA,Ij)", cart_y);
  dpd_buf4_init(&Z1, CC_TMP0, irrep_y, 10, 0, 10, 0, 0, lbl);
  sprintf(lbl, "X_%1s_IjAb (%5.3f)", cart_y, omega_y);
  dpd_buf4_init(&Y2, CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 5, 10, 5, 0, "WAmEf");
  dpd_contract444(&W, &Y2, &Z1, 0, 0, 1, 0);
  dpd_buf4_close(&W);
  dpd_buf4_close(&Y2);

  dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,bA)");
  sprintf(lbl, "X_%1s_IA (%5.3f)", cart_x, omega_x);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);
  dpd_contract244(&X1, &Z1, &Z, 0, 0, 1, -2, 0);
  dpd_file2_close(&X1);
  dpd_buf4_close(&Z1);
  dpd_buf4_sort_axpy(&Z, CC_TMP0, pqsr, 0, 5, "Z(Ij,Ab) Final", 1);
  dpd_buf4_close(&Z);

  dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");
  dpd_buf4_init(&L2, CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
  polar += dpd_buf4_dot(&L2, &Z);
  dpd_buf4_close(&L2);
  dpd_buf4_close(&Z);

  return polar;
}
