#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/* Wmnie_build(): Computes all contributions to the Wmnie HBAR matrix
** elements.  The spin-orbital expression for this term is:
**
** Wmnie = <mn||ie> + t_i^f <mn||fe>
**
** [cf. Gauss and Stanton, JCP 103, 3561-3577 (1995)]
**
** The storage a naming convention for each of the four spin cases
** are as follows:
**
** Spin Case    Storage    Name
** ----------   ---------  -------
** WMNIE        (M>N,EI)   "WMNIE"
** Wmnie        (m>n,ei)   "Wmnie"
** WMnIe        (Mn,eI)    "WMnIe"
** WmNiE        (mN,Ei)    "WmNiE"
** -------------------------------
**
** TDC, June 2002 
*/

void Wmnie_build(void) {
  dpdbuf4 W, Wmnie, WMNIE, WMnIe, WmNiE, WMniE, WmNIe;
  dpdbuf4 E, Z;
  dpdbuf4 D, D_a;
  dpdfile2 t1, tIA, tia;

  if(params.ref == 0) { /** RHF **/

    dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    dpd_buf4_copy(&E, CC_TMP0, "WMnIe (Mn,Ie)");
    dpd_buf4_close(&E);

    /* D(Mn,Fe) * T(I,F) --> W(Mn,Ie) */
    dpd_buf4_init(&WMnIe, CC_TMP0, 0, 0, 10, 0, 10, 0, "WMnIe (Mn,Ie)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract244(&t1, &D, &WMnIe, 1, 2, 1, 1, 1);
    dpd_file2_close(&t1);
    dpd_buf4_close(&D);
    /* W(Mn,Ie) --> W(Mn,eI) */
    dpd_buf4_sort(&WMnIe, CC_HBAR, pqsr, 0, 11, "WMnIe");
    dpd_buf4_close(&WMnIe);

  }
  else if(params.ref == 1) { /** ROHF **/

    /* E(M>N,EI) --> W(M>N,EI) */
    dpd_buf4_init(&E, CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
    dpd_buf4_sort(&E, CC_HBAR, pqsr, 2, 11, "WMNIE");
    dpd_buf4_sort(&E, CC_HBAR, pqsr, 2, 11, "Wmnie");
    dpd_buf4_close(&E);


    /* D(M>N,EF) * T(I,F) --> W(M>N,EI) */
    dpd_buf4_init(&WMNIE, CC_HBAR, 0, 2, 11, 2, 11, 0, "WMNIE");
    dpd_buf4_init(&D_a, CC_DINTS, 0, 2, 5, 2, 5,0, "D <ij||ab> (i>j,ab)");
    dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract424(&D_a,&t1,&WMNIE, 3, 1, 0, -1, 1);
    dpd_file2_close(&t1);
    dpd_buf4_close(&D_a);
    dpd_buf4_close(&WMNIE);


    /* D(m>n,ef) * T(i,f) --> W(m>n,ei) */
    dpd_buf4_init(&Wmnie, CC_HBAR, 0, 2, 11, 2, 11, 0, "Wmnie");
    dpd_buf4_init(&D_a, CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");
    dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tia");
    dpd_contract424(&D_a, &t1, &Wmnie, 3, 1, 0, -1, 1);
    dpd_file2_close(&t1);
    dpd_buf4_close(&D_a);
    dpd_buf4_close(&Wmnie);


    dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    dpd_buf4_copy(&E, CC_TMP0, "WMnIe (Mn,Ie)");
    dpd_buf4_copy(&E, CC_TMP1, "WmNiE (mN,iE)");
    dpd_buf4_close(&E);

    /* D(Mn,Fe) * T(I,F) --> W(Mn,Ie) */
    dpd_buf4_init(&WMnIe, CC_TMP0, 0, 0, 10, 0, 10, 0, "WMnIe (Mn,Ie)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract244(&t1, &D, &WMnIe, 1, 2, 1, 1, 1);
    dpd_file2_close(&t1);
    dpd_buf4_close(&D);
    /* W(Mn,Ie) --> W(Mn,eI) */
    dpd_buf4_sort(&WMnIe, CC_HBAR, pqsr, 0, 11, "WMnIe");
    dpd_buf4_close(&WMnIe);

    /* D(mN,fE) * T(i,f) --> W(mN.iE) */
    dpd_buf4_init(&WmNiE, CC_TMP1, 0, 0, 10, 0, 10, 0, "WmNiE (mN,iE)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tia");
    dpd_contract244(&t1,&D,&WmNiE, 1, 2, 1, 1, 1);
    dpd_file2_close(&t1);
    dpd_buf4_close(&D);
    /* W(mN,iE) --> W(mN,Ei) */
    dpd_buf4_sort(&WmNiE, CC_HBAR, pqsr, 0, 11, "WmNiE");
    dpd_buf4_close(&WmNiE);

    /* also put "normal" sorted versions in CC_HBAR */
    dpd_buf4_init(&WMNIE, CC_HBAR, 0, 2, 11, 2, 11, 0, "WMNIE");
    dpd_buf4_sort(&WMNIE, CC_HBAR, pqsr, 2, 10, "WMNIE (M>N,IE)");
    dpd_buf4_close(&WMNIE);
    dpd_buf4_init(&Wmnie, CC_HBAR, 0, 2, 11, 2, 11, 0, "Wmnie");
    dpd_buf4_sort(&Wmnie, CC_HBAR, pqsr, 2, 10, "Wmnie (m>n,ie)");
    dpd_buf4_close(&Wmnie);
    dpd_buf4_init(&WMnIe, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe");
    dpd_buf4_sort(&WMnIe, CC_HBAR, pqsr, 0, 10, "WMnIe (Mn,Ie)");
    dpd_buf4_close(&WMnIe);
    dpd_buf4_init(&WmNiE, CC_HBAR, 0, 0, 11, 0, 11, 0, "WmNiE");
    dpd_buf4_sort(&WmNiE, CC_HBAR, pqsr, 0, 10, "WmNiE (mN,iE)");
    dpd_buf4_close(&WmNiE);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");

    /* <M>N||IE> --> W(M>N,EI) */
    dpd_buf4_init(&E, CC_EINTS, 0, 2, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
    dpd_buf4_sort(&E, CC_HBAR, pqsr, 2, 21, "WMNIE");
    dpd_buf4_close(&E);

    /* <M>N||EF> T(I,F) --> W(M>N,EI) */
    dpd_buf4_init(&W, CC_HBAR, 0, 2, 21, 2, 21, 0, "WMNIE");
    dpd_buf4_init(&D, CC_DINTS, 0, 2, 5, 2, 5, 0, "D <IJ||AB> (I>J,AB)");
    dpd_contract424(&D, &tIA, &W, 3, 1, 0, -1, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);


    /* <m>n||ie> --> W(m>n,ei) */
    dpd_buf4_init(&E, CC_EINTS, 0, 12, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
    dpd_buf4_sort(&E, CC_HBAR, pqsr, 12, 31, "Wmnie");
    dpd_buf4_close(&E);

    /* <m>n||ef> T(i,f) --> W(m>n,ei) */
    dpd_buf4_init(&W, CC_HBAR, 0, 12, 31, 12, 31, 0, "Wmnie");
    dpd_buf4_init(&D, CC_DINTS, 0, 12, 15, 12, 15, 0, "D <ij||ab> (i>j,ab)");
    dpd_contract424(&D, &tia, &W, 3, 1, 0, -1, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);


    /* <Mn|Ie> --> W(Mn,eI) */
    dpd_buf4_init(&E, CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    dpd_buf4_sort(&E, CC_HBAR, pqsr, 22, 25, "WMnIe");
    dpd_buf4_close(&E);

    /* Z(nM,eI) = <nM|eF> T(I,F) */
    dpd_buf4_init(&Z, CC_TMP1, 0, 23, 25, 23, 25, 0, "Z(nM,eI)");
    dpd_buf4_init(&D, CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    dpd_contract424(&D, &tIA, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&D);
    /* Z(nM,eI) --> W(Mn,eI) */
    dpd_buf4_sort_axpy(&Z, CC_HBAR, qprs, 22, 25, "WMnIe", 1);
    dpd_buf4_close(&Z);


    /* <mN|iE> --> W(mN,Ei) */
    dpd_buf4_init(&E, CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
    dpd_buf4_sort(&E, CC_HBAR, pqsr, 23, 26, "WmNiE");
    dpd_buf4_close(&E);

    /* Z(Nm,Ei) = <Nm|Ef> T(i,f) */
    dpd_buf4_init(&Z, CC_TMP1, 0, 22, 26, 22, 26, 0, "Z(Nm,Ei)");
    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_contract424(&D, &tia, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&D);
    /* Z(Nm,Ei) --> W(mN,Ei) */
    dpd_buf4_sort_axpy(&Z, CC_HBAR, qprs, 23, 26, "WmNiE", 1);
    dpd_buf4_close(&Z);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);

  }

  return;
}
