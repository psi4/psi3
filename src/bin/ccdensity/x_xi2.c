#include <stdio.h>
#include <libdpd/dpd.h>
#include "math.h"
#define EXTERN
#include "globals.h"

extern void x_xi_check(char *term_lbl);
extern void x_xi2_4(void);
extern void x_xi2_14(void);

/* compute xi_2 amplitudes for zeta equations */

void x_xi2(void)
{
  dpdfile2 L1, XIA, Xia, I1, R1, F1, Z1A, Z1B;
  int L_irr, R_irr, G_irr;
  double tval;
  dpdbuf4 D2, R2, L2, H2, I2, Z, Z2, XIJAB, Xijab, XIjAb;

  L_irr = params.L_irr;
  R_irr = params.R_irr;
  G_irr = params.G_irr;

#ifdef DEBUG_XI
  x_xi_check("reset");
#endif

  if(params.ref == 0 || params.ref == 1) {

    /* terms 1 and 5, Xijab += (Lme Rme + 0.25 Lmnef Rmnef) <ij||eb> */
    /* overlaps in params are assigned in x_xi1.c */
    dpd_buf4_init(&D2, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
    dpd_buf4_scmcopy(&D2, EOM_XI, "XIJAB", params.overlap1+params.overlap2);
    dpd_buf4_scmcopy(&D2, EOM_XI, "Xijab", params.overlap1+params.overlap2);
    dpd_buf4_close(&D2);
    dpd_buf4_init(&D2, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_buf4_scmcopy(&D2, EOM_XI, "XIjAb", params.overlap1+params.overlap2);
    dpd_buf4_close(&D2);
#ifdef DEBUG_XI
x_xi_check("terms 1 and 5");
#endif

    /* terms 2 and 9, Xijab -= P(ab) (Lma Rme + Lmnfa Rmnfe) <ij||eb */
    dpd_buf4_init(&Z2, EOM_TMP1, 0, 2, 5, 2, 5, 0, "Z (I>J,AB)");
    dpd_file2_init(&I1, EOM_TMP, G_irr, 1, 1, "LR_VV");
    dpd_buf4_init(&D2, CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");
    dpd_contract244(&I1, &D2, &Z2, 1, 2, 1, 1.0, 0.0);
    dpd_buf4_close(&D2);
    dpd_file2_close(&I1);
    dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 2, 5, 2, 7, 0, "XIJAB");
    dpd_buf4_axpy(&Z2, &XIJAB, -1.0);
    dpd_buf4_close(&XIJAB);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, pqsr, 2, 7, "XIJAB", 1.0);
    dpd_buf4_close(&Z2);

    dpd_buf4_init(&Z2, EOM_TMP1, 0, 2, 5, 2, 5, 0, "Z (i>j,ab)");
    dpd_file2_init(&I1, EOM_TMP, G_irr, 1, 1, "LR_vv");
    dpd_buf4_init(&D2, CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");
    dpd_contract244(&I1, &D2, &Z2, 1, 2, 1, 1.0, 0.0);
    dpd_buf4_close(&D2);
    dpd_file2_close(&I1);
    dpd_buf4_init(&Xijab, EOM_XI, G_irr, 2, 5, 2, 7, 0, "Xijab");
    dpd_buf4_axpy(&Z2, &Xijab, -1.0);
    dpd_buf4_close(&Xijab);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, pqsr, 2, 7, "Xijab", 1.0);
    dpd_buf4_close(&Z2);

    dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
    dpd_file2_init(&I1, EOM_TMP, G_irr, 1, 1, "LR_VV");
    dpd_buf4_init(&D2, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract244(&I1, &D2, &XIjAb, 1, 2, 1, -1.0, 1.0);
    dpd_file2_close(&I1);
    dpd_file2_init(&I1, EOM_TMP, G_irr, 1, 1, "LR_vv");
    dpd_contract424(&D2, &I1, &XIjAb, 3, 1, 0, -1.0, 1.0);
    dpd_file2_close(&I1);
    dpd_buf4_close(&D2);
    dpd_buf4_close(&XIjAb);
#ifdef DEBUG_XI
x_xi_check("terms 2 and 9");
#endif

    /* terms 3 and 10, Xijab -= P(ij) (Lie Rme + 0.5 Linef Rmnef) <mj||ab> */
    dpd_buf4_init(&Z2, EOM_TMP1, 0, 0, 7, 0, 7, 0, "Z (IJ,A>B)");
    dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 0, "LR_OO");
    dpd_buf4_init(&D2, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
    dpd_contract244(&I1, &D2, &Z2, 1, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&D2);
    dpd_file2_close(&I1);
    dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 0, 7, 2, 7, 0, "XIJAB");
    dpd_buf4_axpy(&Z2, &XIJAB, -1.0);
    dpd_buf4_close(&XIJAB);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, qprs, 2, 7, "XIJAB", 1.0);
    dpd_buf4_close(&Z2);

    dpd_buf4_init(&Z2, EOM_TMP1, 0, 0, 7, 0, 7, 0, "Z (ij,a>b)");
    dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 0, "LR_oo");
    dpd_buf4_init(&D2, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
    dpd_contract244(&I1, &D2, &Z2, 1, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&D2);
    dpd_file2_close(&I1);
    dpd_buf4_init(&Xijab, EOM_XI, G_irr, 0, 7, 2, 7, 0, "Xijab");
    dpd_buf4_axpy(&Z2, &Xijab, -1.0);
    dpd_buf4_close(&Xijab);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, qprs, 2, 7, "Xijab", 1.0);
    dpd_buf4_close(&Z2);

    dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
    dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 0, "LR_OO");
    dpd_buf4_init(&D2, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract244(&I1, &D2, &XIjAb, 1, 0, 0, -1.0, 1.0);
    dpd_file2_close(&I1);
    dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 0, "LR_oo");
    dpd_contract424(&D2, &I1, &XIjAb, 1, 1, 1, -1.0, 1.0);
    dpd_file2_close(&I1);
    dpd_buf4_close(&D2);
    dpd_buf4_close(&XIjAb);

    psio_close(EOM_TMP1,0);
    psio_open(EOM_TMP1, PSIO_OPEN_NEW);
#ifdef DEBUG_XI
x_xi_check("terms 3 and 10");
#endif

    x_xi2_4();

    psio_close(EOM_TMP1,0);
    psio_open(EOM_TMP1, PSIO_OPEN_NEW);
#ifdef DEBUG_XI
x_xi_check("terms 4 and 6");
#endif

    /* term 7, Xijab += 0.25 Lmnab Rmnef <ij||ef> */
    dpd_buf4_init(&D2, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
    dpd_buf4_init(&R2, CC_GR, R_irr, 2, 7, 2, 7, 0, "RIJAB");
    dpd_buf4_init(&Z2, EOM_TMP1, R_irr, 2, 2, 2, 2, 0, "Z (I>J,M>N)");
    dpd_contract444(&D2, &R2, &Z2, 0, 0, 1.0, 0.0); 
    dpd_buf4_close(&R2);
    dpd_buf4_close(&D2);
    dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 2, 7, 2, 7, 0, "XIJAB");
    dpd_buf4_init(&L2, CC_GL, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    dpd_contract444(&Z2, &L2, &XIJAB, 0, 1, 1.0, 1.0); 
    dpd_buf4_close(&L2);
    dpd_buf4_close(&Z2);
    dpd_buf4_close(&XIJAB);

    dpd_buf4_init(&D2, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
    dpd_buf4_init(&R2, CC_GR, R_irr, 2, 7, 2, 7, 0, "Rijab");
    dpd_buf4_init(&Z2, EOM_TMP1, R_irr, 2, 2, 2, 2, 0, "Z (i>j,m>n)");
    dpd_contract444(&D2, &R2, &Z2, 0, 0, 1.0, 0.0); 
    dpd_buf4_close(&R2);
    dpd_buf4_close(&D2);
    dpd_buf4_init(&Xijab, EOM_XI, G_irr, 2, 7, 2, 7, 0, "Xijab");
    dpd_buf4_init(&L2, CC_GL, L_irr, 2, 7, 2, 7, 0, "Lijab");
    dpd_contract444(&Z2, &L2, &Xijab, 0, 1, 1.0, 1.0); 
    dpd_buf4_close(&L2);
    dpd_buf4_close(&Z2);
    dpd_buf4_close(&Xijab);

    dpd_buf4_init(&D2, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_buf4_init(&R2, CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjAb");
    dpd_buf4_init(&Z2, EOM_TMP1, R_irr, 0, 0, 0, 0, 0, "Z (Ij,Mn)");
    dpd_contract444(&D2, &R2, &Z2, 0, 0, 1.0, 0.0); 
    dpd_buf4_close(&R2);
    dpd_buf4_close(&D2);
    dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
    dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_contract444(&Z2, &L2, &XIjAb, 0, 1, 1.0, 1.0); 
    dpd_buf4_close(&L2);
    dpd_buf4_close(&Z2);
    dpd_buf4_close(&XIjAb);

    psio_close(EOM_TMP1,0);
    psio_open(EOM_TMP1, PSIO_OPEN_NEW);
#ifdef DEBUG_XI
x_xi_check("term 7");
#endif

    /* term 8, Xijab += 0.25 Rmnef Lijef <mn||ab> */
    dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 2, 7, 2, 7, 0, "XIJAB");
    dpd_buf4_init(&I2, EOM_TMP, G_irr, 2, 2, 2, 2, 0, "R2L2_OOOO");
    dpd_buf4_init(&D2, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
    dpd_contract444(&I2, &D2, &XIJAB, 1, 1, 1.0, 1.0); 
    dpd_buf4_close(&D2);
    dpd_buf4_close(&I2);
    dpd_buf4_close(&XIJAB);
    dpd_buf4_init(&Xijab, EOM_XI, G_irr, 2, 7, 2, 7, 0, "Xijab");
    dpd_buf4_init(&I2, EOM_TMP, G_irr, 2, 2, 2, 2, 0, "R2L2_oooo");
    dpd_buf4_init(&D2, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
    dpd_contract444(&I2, &D2, &Xijab, 1, 1, 1.0, 1.0); 
    dpd_buf4_close(&D2);
    dpd_buf4_close(&I2);
    dpd_buf4_close(&Xijab);
    dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
    dpd_buf4_init(&I2, EOM_TMP, G_irr, 0, 0, 0, 0, 0, "R2L2_OoOo");
    dpd_buf4_init(&D2, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract444(&I2, &D2, &XIjAb, 1, 1, 1.0, 1.0); 
    dpd_buf4_close(&D2);
    dpd_buf4_close(&I2);
    dpd_buf4_close(&XIjAb);

#ifdef DEBUG_XI
x_xi_check("term 8");
#endif

    /* term 11, Xijab -= 0.5 P(ab) Lijfb (Rmnef <mn||ea>) */
    /* term 17        -=     P(ab) Lijfb (Rmf Fma) */
    /* term 20        +=     P(ab) Lijfb (Rme Wfmae) */
    /* build 1-e intermediates to include term 17 */
        /* for term 11: */
    dpd_file2_init(&I1, EOM_TMP_XI, R_irr, 1, 1, "RD_VV");
    dpd_file2_copy(&I1, EOM_TMP1, "Z (F,A)");
    dpd_file2_close(&I1);
    dpd_file2_init(&I1, EOM_TMP1, R_irr, 1, 1, "Z (F,A)");
        /* for term 20: */
    dpd_file2_init(&Z1A, EOM_TMP_XI, R_irr, 1, 1, "R1Wamef_VV");
    dpd_file2_axpy(&Z1A, &I1, -1.0, 0);
    dpd_file2_close(&Z1A);
        /* for term 17: */
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_file2_init(&F1, CC_OEI, 0, 0, 1, "FME");
    dpd_contract222(&R1, &F1, &I1, 1, 1, 1.0, 1.0);
    dpd_file2_close(&F1);
    dpd_file2_close(&R1);
    dpd_file2_close(&I1);

        /* for term 11: */
    dpd_file2_init(&I1, EOM_TMP_XI, R_irr, 1, 1, "RD_vv");
    dpd_file2_copy(&I1, EOM_TMP1, "Z (f,a)");
    dpd_file2_close(&I1);
    dpd_file2_init(&I1, EOM_TMP1, R_irr, 1, 1, "Z (f,a)");
        /* for term 20: */
    dpd_file2_init(&Z1B, EOM_TMP_XI, R_irr, 1, 1, "R1Wamef_vv");
    dpd_file2_axpy(&Z1B, &I1, -1.0, 0);
    dpd_file2_close(&Z1B);
        /* for term 17: */
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
    dpd_file2_init(&F1, CC_OEI, 0, 0, 1, "Fme");
    dpd_contract222(&R1, &F1, &I1, 1, 1, 1.0, 1.0);
    dpd_file2_close(&F1);
    dpd_file2_close(&R1);
    dpd_file2_close(&I1);

    dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 2, 5, 2, 5, 0, "Z (I>J,AB)"); 
    dpd_file2_init(&I1, EOM_TMP1, R_irr, 1, 1, "Z (F,A)");
    dpd_buf4_init(&L2, CC_GL, L_irr, 2, 5, 2, 7, 0, "LIJAB");
    dpd_contract244(&I1, &L2, &Z2, 0, 2, 1, 1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_file2_close(&I1);
    dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 2, 5, 2, 7, 0, "XIJAB");
    dpd_buf4_axpy(&Z2, &XIJAB, -1.0);
    dpd_buf4_close(&XIJAB);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, pqsr, 2, 7, "XIJAB", 1.0);
    dpd_buf4_close(&Z2);

    dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 2, 5, 2, 5, 0, "Z (i>j,ab)"); 
    dpd_file2_init(&I1, EOM_TMP1, R_irr, 1, 1, "Z (f,a)");
    dpd_buf4_init(&L2, CC_GL, L_irr, 2, 5, 2, 7, 0, "Lijab");
    dpd_contract244(&I1, &L2, &Z2, 0, 2, 1, 1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_file2_close(&I1);
    dpd_buf4_init(&Xijab, EOM_XI, G_irr, 2, 5, 2, 7, 0, "Xijab");
    dpd_buf4_axpy(&Z2, &Xijab, -1.0);
    dpd_buf4_close(&Xijab);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, pqsr, 2, 7, "Xijab", 1.0);
    dpd_buf4_close(&Z2);

    dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
    dpd_file2_init(&I1, EOM_TMP1, R_irr, 1, 1, "Z (F,A)");
    dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_contract244(&I1, &L2, &XIjAb, 0, 2, 1, -1.0, 1.0);
    dpd_file2_close(&I1);
    dpd_file2_init(&I1, EOM_TMP1, R_irr, 1, 1, "Z (f,a)");
    dpd_contract424(&L2, &I1, &XIjAb, 3, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_file2_close(&I1);
    dpd_buf4_close(&XIjAb);
#ifdef DEBUG_XI
x_xi_check("terms 11, 17 and 20");
#endif

    /* term 12, Xijab -= 0.5 P(ij) Lmjab (Rmnef <in||ef>) */
    /* term 16,       -=     P(ij) Lmjab (Rme Fie) */
    /* term 21,       -=     P(ij) Lmjab (Rne Winme) */
    /* make 1-electron intermediates to include terms 16 and 21 as well */
        /* for term 12: */
    dpd_file2_init(&I1, EOM_TMP_XI, R_irr, 0, 0, "RD_OO");
    dpd_file2_copy(&I1, EOM_TMP1, "Z (M,I)");
    dpd_file2_close(&I1);
    dpd_file2_init(&I1, EOM_TMP1, R_irr, 0, 0, "Z (M,I)");
         /* for term 21 */
    dpd_file2_init(&Z1A, EOM_TMP_XI, R_irr, 0, 0, "R1Wmnie_OO");
    dpd_file2_axpy(&Z1A, &I1, 1.0, 1);
    dpd_file2_close(&Z1A);
        /* for term 16: */
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_file2_init(&F1, CC_OEI, 0, 0, 1, "FME");
    dpd_contract222(&R1, &F1, &I1, 0, 0, 1.0, 1.0);
    dpd_file2_close(&F1);
    dpd_file2_close(&R1);
    dpd_file2_close(&I1);

        /* for term 12 */
    dpd_file2_init(&I1, EOM_TMP_XI, R_irr, 0, 0, "RD_oo");
    dpd_file2_copy(&I1, EOM_TMP1, "Z (m,i)");
    dpd_file2_close(&I1);
    dpd_file2_init(&I1, EOM_TMP1, R_irr, 0, 0, "Z (m,i)");
        /* for term 21 */
    dpd_file2_init(&Z1B, EOM_TMP_XI, R_irr, 0, 0, "R1Wmnie_oo");
    dpd_file2_axpy(&Z1B, &I1, 1.0, 1);
    dpd_file2_close(&Z1B);
        /* for term 16 */
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
    dpd_file2_init(&F1, CC_OEI, 0, 0, 1, "Fme");
    dpd_contract222(&R1, &F1, &I1, 0, 0, 1.0, 1.0);
    dpd_file2_close(&F1);
    dpd_file2_close(&R1);
    dpd_file2_close(&I1);

    dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 0, 7, 0, 7, 0, "Z (IJ,A>B)"); 
    dpd_file2_init(&I1, EOM_TMP1, R_irr, 0, 0, "Z (M,I)");
    dpd_buf4_init(&L2, CC_GL, L_irr, 0, 7, 2, 7, 0, "LIJAB");
    dpd_contract244(&I1, &L2, &Z2, 0, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_file2_close(&I1);
    dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 0, 7, 2, 7, 0, "XIJAB");
    dpd_buf4_axpy(&Z2, &XIJAB, -1.0);
    dpd_buf4_close(&XIJAB);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, qprs, 2, 7, "XIJAB", 1.0);
    dpd_buf4_close(&Z2);

    dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 0, 7, 0, 7, 0, "Z (ij,a>b)"); 
    dpd_file2_init(&I1, EOM_TMP1, R_irr, 0, 0, "Z (m,i)");
    dpd_buf4_init(&L2, CC_GL, L_irr, 0, 7, 2, 7, 0, "Lijab");
    dpd_contract244(&I1, &L2, &Z2, 0, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_file2_close(&I1);
    dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 0, 7, 2, 7, 0, "Xijab");
    dpd_buf4_axpy(&Z2, &XIJAB, -1.0);
    dpd_buf4_close(&XIJAB);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, qprs, 2, 7, "Xijab", 1.0);
    dpd_buf4_close(&Z2);

    dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
    dpd_file2_init(&I1, EOM_TMP1, R_irr, 0, 0, "Z (M,I)");
    dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_contract244(&I1, &L2, &XIjAb, 0, 0, 0, -1.0, 1.0);
    dpd_file2_close(&I1);
    dpd_file2_init(&I1, EOM_TMP1, R_irr, 0, 0, "Z (m,i)");
    dpd_contract424(&L2, &I1, &XIjAb, 1, 0, 1, -1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_file2_close(&I1);
    dpd_buf4_close(&XIjAb);
#ifdef DEBUG_XI
x_xi_check("term 12, 16 and 21");
#endif

    psio_close(EOM_TMP1,0);
    psio_open(EOM_TMP1, PSIO_OPEN_NEW);

    /* term 13 + 15, (0.25 Rmnef <mn||ef> + Rme Fme) Lijab */
    if (L_irr == 0) {
      dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
      dpd_file2_init(&F1, CC_OEI, 0, 0, 1, "FME");
      tval = dpd_file2_dot(&R1, &F1);
      dpd_file2_close(&F1);
      dpd_file2_close(&R1);
      dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
      dpd_file2_init(&F1, CC_OEI, 0, 0, 1, "Fme");
      tval += dpd_file2_dot(&R1, &F1);
      dpd_file2_close(&F1);
      dpd_file2_close(&R1);
      tval += params.RD_overlap;
  
      dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 2, 7, 2, 7, 0, "XIJAB");
      dpd_buf4_init(&L2, CC_GL, L_irr, 2, 7, 2, 7, 0, "LIJAB");
      dpd_buf4_axpy(&L2, &XIJAB, tval);
      dpd_buf4_close(&L2);
      dpd_buf4_close(&XIJAB);
      dpd_buf4_init(&Xijab, EOM_XI, G_irr, 2, 7, 2, 7, 0, "Xijab");
      dpd_buf4_init(&L2, CC_GL, L_irr, 2, 7, 2, 7, 0, "Lijab");
      dpd_buf4_axpy(&L2, &Xijab, tval);
      dpd_buf4_close(&L2);
      dpd_buf4_close(&Xijab);
      dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
      dpd_buf4_init(&L2, CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
      dpd_buf4_axpy(&L2, &XIjAb, tval);
      dpd_buf4_close(&L2);
      dpd_buf4_close(&XIjAb);
#ifdef DEBUG_XI
x_xi_check("term 13 (ijab) and 15 (Fme)");
#endif
    }

    /* term 14, +P(ij) P(ab) Lmjeb Rme Fia */
    x_xi2_14();
#ifdef DEBUG_XI
x_xi_check("term 14 (Fme)");
#endif
    psio_close(EOM_TMP1,0);
    psio_open(EOM_TMP1, PSIO_OPEN_NEW);

    /* term 22, +P(ij) (Limfe Rme) Wfjab */
    dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 0, 7, 0, 7, 0, "Z (IJ,A>B)");
    dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    dpd_buf4_init(&H2, CC_HBAR, 0, 11, 7, 11, 7, 0, "WAMEF (AM,E>F)");
    dpd_contract244(&I1, &H2, &Z2, 1, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&H2);
    dpd_file2_close(&I1);
    dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 0, 7, 2, 7, 0, "XIJAB");
    dpd_buf4_axpy(&Z2, &XIJAB, 1.0);
    dpd_buf4_close(&XIJAB);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, qprs, 2, 7, "XIJAB", -1.0);
    dpd_buf4_close(&Z2);

    dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 0, 7, 0, 7, 0, "Z (ij,a>b)");
    dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_ov");
    dpd_buf4_init(&H2, CC_HBAR, 0, 11, 7, 11, 7, 0, "Wamef (am,e>f)");
    dpd_contract244(&I1, &H2, &Z2, 1, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&H2);
    dpd_file2_close(&I1);
    dpd_buf4_init(&Xijab, EOM_XI, G_irr, 0, 7, 2, 7, 0, "Xijab");
    dpd_buf4_axpy(&Z2, &Xijab, 1.0);
    dpd_buf4_close(&Xijab);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, qprs, 2, 7, "Xijab", -1.0);
    dpd_buf4_close(&Z2);

    dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
    dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf (Am,Ef)");
    dpd_contract244(&I1, &H2, &XIjAb, 1, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&H2);
    dpd_file2_close(&I1);
    dpd_buf4_close(&XIjAb);

    dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z (jI,bA)");
    dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_ov");
    dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 5, 0, "WaMeF (aM,eF)");
    dpd_contract244(&I1, &H2, &Z2, 1, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&H2);
    dpd_file2_close(&I1);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, qpsr, 0, 5, "XIjAb", 1.0);
    dpd_buf4_close(&Z2);
#ifdef DEBUG_XI
x_xi_check("term 22 (Wamef)");
#endif
    psio_close(EOM_TMP1,0);
    psio_open(EOM_TMP1, PSIO_OPEN_NEW);

    /* term 23, -P(ab) (Lmnea Rme) Wijnb */
    dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 2, 5, 2, 5, 0, "Z (I>J,AB)");
    dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    dpd_buf4_init(&H2, CC_HBAR, 0, 2, 10, 2, 10, 0, "WMNIE (M>N,IE)");
    dpd_contract244(&I1, &H2, &Z2, 0, 2, 1, 1.0, 0.0);
    dpd_buf4_close(&H2);
    dpd_file2_close(&I1);
    dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 2, 5, 2, 7, 0, "XIJAB");
    dpd_buf4_axpy(&Z2, &XIJAB, -1.0);
    dpd_buf4_close(&XIJAB);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, pqsr, 2, 7, "XIJAB", 1.0);
    dpd_buf4_close(&Z2);

    dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 2, 5, 2, 5, 0, "Z (i>j,ab)");
    dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_ov");
    dpd_buf4_init(&H2, CC_HBAR, 0, 2, 10, 2, 10, 0, "Wmnie (m>n,ie)");
    dpd_contract244(&I1, &H2, &Z2, 0, 2, 1, 1.0, 0.0);
    dpd_buf4_close(&H2);
    dpd_file2_close(&I1);
    dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 2, 5, 2, 7, 0, "Xijab");
    dpd_buf4_axpy(&Z2, &XIJAB, -1.0);
    dpd_buf4_close(&XIJAB);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, pqsr, 2, 7, "Xijab", 1.0);
    dpd_buf4_close(&Z2);

    dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
    dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    dpd_buf4_init(&H2, CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe (Mn,Ie)");
    dpd_contract244(&I1, &H2, &XIjAb, 0, 2, 1, -1.0, 1.0);
    dpd_buf4_close(&H2);
    dpd_file2_close(&I1);
    dpd_buf4_close(&XIjAb);

    dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z (jI,Ab)");
    dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_ov");
    dpd_buf4_init(&H2, CC_HBAR, 0, 0, 11, 0, 11, 0, "WmNiE");
    dpd_contract424(&H2, &I1, &Z2, 3, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&H2);
    dpd_file2_close(&I1);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, qprs, 0, 5, "XIjAb", -1.0);
    dpd_buf4_close(&Z2);
#ifdef DEBUG_XI
x_xi_check("term 23 (Wmnie)");
#endif

    /* term 25, Xijab += (Lnmab Rme) Wijne */
    dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 2, 7, 2, 7, 0, "XIJAB");
    dpd_buf4_init(&I2, EOM_TMP, G_irr, 7, 10, 7, 10, 0, "L2R1_VVOV");
    dpd_buf4_init(&H2, CC_HBAR, 0, 2, 10, 2, 10, 0, "WMNIE (M>N,IE)");
    dpd_contract444(&H2, &I2, &XIJAB, 0, 0, 1.0, 1.0); 
    dpd_buf4_close(&H2);
    dpd_buf4_close(&I2);
    dpd_buf4_close(&XIJAB);

    dpd_buf4_init(&Xijab, EOM_XI, G_irr, 2, 7, 2, 7, 0, "Xijab");
    dpd_buf4_init(&I2, EOM_TMP, G_irr, 7, 10, 7, 10, 0, "L2R1_vvov");
    dpd_buf4_init(&H2, CC_HBAR, 0, 2, 10, 2, 10, 0, "Wmnie (m>n,ie)");
    dpd_contract444(&H2, &I2, &Xijab, 0, 0, 1.0, 1.0); 
    dpd_buf4_close(&H2);
    dpd_buf4_close(&I2);
    dpd_buf4_close(&Xijab);

    dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
    dpd_buf4_init(&I2, EOM_TMP, G_irr, 5, 10, 5, 10, 0, "L2R1_VvOv");
    dpd_buf4_init(&H2, CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe (Mn,Ie)");
    dpd_contract444(&H2, &I2, &XIjAb, 0, 0, 1.0, 1.0); 
    dpd_buf4_close(&H2);
    dpd_buf4_close(&I2);
    dpd_buf4_close(&XIjAb);

    dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z (jI,Ab)");
    dpd_buf4_init(&I2, EOM_TMP, G_irr, 5, 10, 5, 10, 0, "L2R1_VvoV");
    dpd_buf4_init(&H2, CC_HBAR, 0, 0, 10, 0, 10, 0, "WmNiE (mN,iE)");
    dpd_contract444(&H2, &I2, &Z2, 0, 0, 1.0, 0.0); 
    dpd_buf4_close(&H2);
    dpd_buf4_close(&I2);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, qprs, 0, 5, "XIjAb", 1.0);
    dpd_buf4_close(&Z2);

#ifdef DEBUG_XI
x_xi_check("term 25 (Wmnie)");
#endif
    psio_close(EOM_TMP1,0);
    psio_open(EOM_TMP1, PSIO_OPEN_NEW);

    /* term 24, Xijab -= (Lijfe Rme) Wfmab */
    dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 2, 7, 2, 7, 0, "XIJAB");
    dpd_buf4_init(&I2, EOM_TMP, G_irr, 2, 11, 2, 11, 0, "L2R1_OOVO");
    dpd_buf4_init(&H2, CC_HBAR, 0, 11, 7, 11, 7, 0, "WAMEF (AM,E>F)");
    dpd_contract444(&I2, &H2, &XIJAB, 0, 1, -1.0, 1.0); 
    dpd_buf4_close(&H2);
    dpd_buf4_close(&I2);
    dpd_buf4_close(&XIJAB);

    dpd_buf4_init(&Xijab, EOM_XI, G_irr, 2, 7, 2, 7, 0, "Xijab");
    dpd_buf4_init(&I2, EOM_TMP, G_irr, 2, 11, 2, 11, 0, "L2R1_oovo");
    dpd_buf4_init(&H2, CC_HBAR, 0, 11, 7, 11, 7, 0, "Wamef (am,e>f)");
    dpd_contract444(&I2, &H2, &Xijab, 0, 1, -1.0, 1.0); 
    dpd_buf4_close(&H2);
    dpd_buf4_close(&I2);
    dpd_buf4_close(&Xijab);

    dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
    dpd_buf4_init(&I2, EOM_TMP, G_irr, 0, 11, 0, 11, 0, "L2R1_OoVo");
    dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf (Am,Ef)");
    dpd_contract444(&I2, &H2, &XIjAb, 0, 1, -1.0, 1.0); 
    dpd_buf4_close(&H2);
    dpd_buf4_close(&I2);
    dpd_buf4_close(&XIjAb);

    dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z (Ij,bA)");
    dpd_buf4_init(&I2, EOM_TMP, G_irr, 0, 11, 0, 11, 0, "L2R1_OovO");
    dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 5, 0, "WaMeF (aM,eF)");
    dpd_contract444(&I2, &H2, &Z2, 0, 1, 1.0, 0.0); 
    dpd_buf4_close(&H2);
    dpd_buf4_close(&I2);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, pqsr, 0, 5, "XIjAb", -1.0);
    dpd_buf4_close(&Z2);
#ifdef DEBUG_XI
x_xi_check("term 24 (Wamef)");
#endif
    psio_close(EOM_TMP1,0);
    psio_open(EOM_TMP1, PSIO_OPEN_NEW);

    /* terms 18, 19: Xijab -= P(ij) P(ab) Linae (Rme Wmjnb + Rnf Wejbf) */
    /* construct Z(JB,NE) = RME WMJNB + RNF WEJBF */
    dpd_buf4_init(&Z, EOM_TMP1, R_irr, 11, 11, 11, 11, 0, "Z (EJ,BN)");
    dpd_buf4_init(&H2, CC_HBAR, 0, 0, 11, 2, 11, 0, "WMNIE");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_contract244(&R1, &H2, &Z, 0, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&H2);
    dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 7, 0, "WAMEF (AM,E>F)");
    dpd_contract424(&H2, &R1, &Z, 3, 1, 0, 1.0, 1.0);
    dpd_buf4_close(&H2);
    dpd_file2_close(&R1);
    dpd_buf4_sort(&Z, EOM_TMP1, qpsr, 10, 10, "Z (JE,NB)");
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (JE,NB)");
    dpd_buf4_sort(&Z, EOM_TMP1, psrq, 10, 10, "Z (JB,NE)");
    dpd_buf4_close(&Z);

    dpd_buf4_init(&Z, EOM_TMP1, R_irr, 11, 10, 11, 10, 0, "Z (eJ,nB)");
    dpd_buf4_init(&H2, CC_HBAR, 0, 0, 10, 0, 10, 0, "WmNiE (mN,iE)");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
    dpd_contract244(&R1, &H2, &Z, 0, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&H2);
    dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 5, 0, "WaMeF (aM,eF)");
    dpd_contract244(&R1, &H2, &Z, 1, 2, 1, -1.0, 1.0);
    dpd_buf4_close(&H2);
    dpd_file2_close(&R1);
    dpd_buf4_sort(&Z, EOM_TMP1, qprs, 10, 10, "Z (Je,nB)");
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (Je,nB)");
    dpd_buf4_sort(&Z, EOM_TMP1, psrq, 10, 10, "Z (JB,ne)");
    dpd_buf4_close(&Z);

    dpd_buf4_init(&Z, EOM_TMP1, R_irr, 11, 11, 11, 11, 0, "Z (ej,bn)");
    dpd_buf4_init(&H2, CC_HBAR, 0, 0, 11, 2, 11, 0, "Wmnie");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
    dpd_contract244(&R1, &H2, &Z, 0, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&H2);
    dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 7, 0, "Wamef (am,e>f)");
    dpd_contract424(&H2, &R1, &Z, 3, 1, 0, 1.0, 1.0);
    dpd_buf4_close(&H2);
    dpd_file2_close(&R1);
    dpd_buf4_sort(&Z, EOM_TMP1, qpsr, 10, 10, "Z (je,nb)");
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (je,nb)");
    dpd_buf4_sort(&Z, EOM_TMP1, psrq, 10, 10, "Z (jb,ne)");
    dpd_buf4_close(&Z);

    dpd_buf4_init(&Z, EOM_TMP1, R_irr, 11, 10, 11, 10, 0, "Z (Ej,Nb)");
    dpd_buf4_init(&H2, CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe (Mn,Ie)");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_contract244(&R1, &H2, &Z, 0, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&H2);
    dpd_buf4_init(&H2, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf (Am,Ef)");
    dpd_contract244(&R1, &H2, &Z, 1, 2, 1, -1.0, 1.0);
    dpd_buf4_close(&H2);
    dpd_file2_close(&R1);
    dpd_buf4_sort(&Z, EOM_TMP1, qprs, 10, 10, "Z (jE,Nb)");
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (jE,Nb)");
    dpd_buf4_sort(&Z, EOM_TMP1, psrq, 10, 10, "Z (jb,NE)");
    dpd_buf4_close(&Z);

    dpd_buf4_init(&Z, EOM_TMP1, R_irr, 10, 11, 10, 11, 0, "Z (Je,bN)");
    dpd_buf4_init(&H2, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
    dpd_contract424(&H2, &R1, &Z, 1, 0, 1, -1.0, 0.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&H2);
    dpd_buf4_init(&H2, CC_HBAR, 0, 10, 5, 10, 5, 0, "WaMeF");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_contract424(&H2, &R1, &Z, 3, 1, 0, 1.0, 1.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&H2);
    dpd_buf4_sort(&Z, EOM_TMP1, prsq, 10, 10, "Z (Jb,Ne)");
    dpd_buf4_close(&Z);

    dpd_buf4_init(&Z, EOM_TMP1, R_irr, 10, 11, 10, 11, 0, "Z (jE,Bn)");
    dpd_buf4_init(&H2, CC_HBAR, 0, 0, 11, 0, 11, 0, "WmNiE");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_contract424(&H2, &R1, &Z, 1, 0, 1, -1.0, 0.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&H2);
    dpd_buf4_init(&H2, CC_HBAR, 0, 10, 5, 10, 5, 0, "WAmEf");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
    dpd_contract424(&H2, &R1, &Z, 3, 1, 0, 1.0, 1.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&H2);
    dpd_buf4_sort(&Z, EOM_TMP1, prsq, 10, 10, "Z (jB,nE)");
    dpd_buf4_close(&Z);

    /* XIJAB -= P(IJ) P(AB) L(IA,NE) Z(NE,JB) */
    dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z2 (IA,JB)");
    dpd_buf4_init(&L2, CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAJB");
    dpd_buf4_init(&Z, EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (JB,NE)");
    dpd_contract444(&L2, &Z, &Z2, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAjb");
    dpd_buf4_init(&Z, EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (JB,ne)");
    dpd_contract444(&L2, &Z, &Z2, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&L2);
    dpd_buf4_sort(&Z2, EOM_TMP1, prqs, 0, 5, "Z2 (IJ,AB)");
    dpd_buf4_close(&Z2);

    dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z2 (IJ,AB)");
    dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 0, 5, 2, 7, 0, "XIJAB");
    dpd_buf4_axpy(&Z2, &XIJAB, -1.0);
    dpd_buf4_close(&XIJAB);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, qprs, 2, 7, "XIJAB", 1.0);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, pqsr, 2, 7, "XIJAB", 1.0);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, qpsr, 2, 7, "XIJAB", -1.0);
    dpd_buf4_close(&Z2);

#ifdef DEBUG_XI
x_xi_check("terms 18, 19 (Wmnie, Wamef)");
#endif

    /* Xijab -= P(ij) P(ab) L(ia,ne) Z(ne,jb) */
    dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z2 (ia,jb)");
    dpd_buf4_init(&L2, CC_GL, L_irr, 10, 10, 10, 10, 0, "Liajb");
    dpd_buf4_init(&Z, EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (jb,ne)");
    dpd_contract444(&L2, &Z, &Z2, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GL, L_irr, 10, 10, 10, 10, 0, "LiaJB");
    dpd_buf4_init(&Z, EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (jb,NE)");
    dpd_contract444(&L2, &Z, &Z2, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&L2);
    dpd_buf4_sort(&Z2, EOM_TMP1, prqs, 0, 5, "Z2 (ij,ab)");
    dpd_buf4_close(&Z2);

    dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z2 (ij,ab)");
    dpd_buf4_init(&Xijab, EOM_XI, G_irr, 0, 5, 2, 7, 0, "Xijab");
    dpd_buf4_axpy(&Z2, &Xijab, -1.0);
    dpd_buf4_close(&Xijab);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, qprs, 2, 7, "Xijab", 1.0);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, pqsr, 2, 7, "Xijab", 1.0);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, qpsr, 2, 7, "Xijab", -1.0);
    dpd_buf4_close(&Z2);

    /* XIjAb += - L(IA,NE) Z(jb,NE) - L(IA,ne) Z(jb,ne)  */
    /*          - L(jb,ne) Z(IA,ne) - L(jb,NE) Z(IA,NE)  */
    dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z2 (IA,jb)");
    dpd_buf4_init(&L2, CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAJB");
    dpd_buf4_init(&Z, EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (jb,NE)");
    dpd_contract444(&L2, &Z, &Z2, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAjb");
    dpd_buf4_init(&Z, EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (jb,ne)");
    dpd_contract444(&L2, &Z, &Z2, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&Z, EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (JB,ne)");
    dpd_buf4_init(&L2, CC_GL, L_irr, 10, 10, 10, 10, 0, "Liajb");
    dpd_contract444(&Z, &L2, &Z2, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (JB,NE)");
    dpd_buf4_init(&L2, CC_GL, L_irr, 10, 10, 10, 10, 0, "LiaJB");
    dpd_contract444(&Z, &L2, &Z2, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&Z);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, prqs, 0, 5, "XIjAb", -1.0);
    dpd_buf4_close(&Z2);

    /* XIjAb += + L(jA,Ne) Z(Ib,Ne) + L(Ib,nE) Z(jA,nE)  */
    dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z2 (Ib,jA)");
    dpd_buf4_init(&L2, CC_GL, L_irr, 10, 10, 10, 10, 0, "LjAIb");
    dpd_buf4_init(&Z, EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (Jb,Ne)");
    dpd_contract444(&Z, &L2, &Z2, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_GL, L_irr, 10, 10, 10, 10, 0, "LIbjA");
    dpd_buf4_init(&Z, EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (jB,nE)");
    dpd_contract444(&L2, &Z, &Z2, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&L2);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, prsq, 0, 5, "XIjAb", 1.0);
    dpd_buf4_close(&Z2);
#ifdef DEBUG_XI
x_xi_check("terms 18, 19 (Wmnie, Wamef)");
#endif
    psio_close(EOM_TMP1,0);
    psio_open(EOM_TMP1, PSIO_OPEN_NEW);
  }

  /* Write irrep of XI amplitudes to CC_INFO */
  psio_write_entry(CC_INFO, "XI Irrep", (char *) &G_irr,sizeof(int));

  dpd_file2_init(&XIA, EOM_XI, G_irr, 0, 1, "XIA");
  tval = dpd_file2_dot_self(&XIA);
  dpd_file2_close(&XIA);
  fprintf(outfile,"XIA amplitudes: norm=%20.15lf\n", sqrt(tval) );
  dpd_file2_init(&Xia, EOM_XI, G_irr, 0, 1, "Xia");
  tval += dpd_file2_dot_self(&Xia);
  dpd_file2_close(&Xia);
  dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 2, 7, 2, 7, 0, "XIJAB");
  tval += dpd_buf4_dot_self(&XIJAB);
  dpd_buf4_close(&XIJAB);
  dpd_buf4_init(&Xijab, EOM_XI, G_irr, 2, 7, 2, 7, 0, "Xijab");
  tval += dpd_buf4_dot_self(&Xijab);
  dpd_buf4_close(&Xijab);
  dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  tval += dpd_buf4_dot_self(&XIjAb);
  dpd_buf4_close(&XIjAb);
  fprintf(outfile,"Norm of Xi: %20.15lf\n", sqrt(tval) );
  return;
}
