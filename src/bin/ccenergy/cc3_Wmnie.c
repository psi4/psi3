#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/* cc3_Wmnie(): Compute the Wmnie matrix from CC3 theory, which is
** given in spin-orbitals as:
** 
** Wmnie = <mn||ie> + t_i^f <mn||fe>
**
** TDC, Feb 2004
*/

extern void purge_HET1_Wmnie(void);

void cc3_Wmnie(void)
{
  dpdbuf4 E, D, W, Z;
  dpdfile2 tIA, tia;

  if(params.ref == 0) { /** RHF **/

    dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    dpd_buf4_copy(&E, CC3_HET1, "CC3 WMnIe (Mn,Ie)");
    dpd_buf4_close(&E);

    dpd_buf4_init(&W, CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WMnIe (Mn,Ie)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract244(&tIA, &D, &W, 1, 2, 1, 1, 1);
    dpd_file2_close(&tIA);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);
  }

  else if (params.ref == 1) { /* ROHF */

    /** W(M>N,IE) <--- <MN||IE> **/
    /** W(m>n,ie) <--- <mn||ie> **/
    dpd_buf4_init(&E, CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
    dpd_buf4_sort(&E, CC3_HET1, pqsr, 2, 11, "CC3 WMNIE (M>N,EI)");
    dpd_buf4_sort(&E, CC3_HET1, pqsr, 2, 11, "CC3 Wmnie (m>n,ei)");
    dpd_buf4_close(&E);

    /** W(Mn,Ie) <--- <Mn|Ie> **/
    /** W(mN,iE) <--- <mN|iE> **/
    dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    dpd_buf4_sort(&E, CC3_HET1, pqsr, 0, 11, "CC3 WMnIe (Mn,eI)");
    dpd_buf4_sort(&E, CC3_HET1, pqsr, 0, 11, "CC3 WmNiE (mN,Ei)");
    dpd_buf4_close(&E);

    /**** Term 2 ****/

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    /* <M>N||EF> T(I,F) --> W(M>N,EI) */
    /* <m>n||ef> T(i,f) --> W(m>n,ei) */
    dpd_buf4_init(&W, CC3_HET1, 0, 2, 11, 2, 11, 0, "CC3 WMNIE (M>N,EI)");
    dpd_buf4_init(&D, CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");
    dpd_contract424(&D, &tIA, &W, 3, 1, 0, -1, 1.0);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC3_HET1, 0, 2, 11, 2, 11, 0, "CC3 Wmnie (m>n,ei)");
    dpd_contract424(&D, &tia, &W, 3, 1, 0, -1, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);

    /* Z(nM,eI) = <nM|eF> T(I,F) */
    dpd_buf4_init(&Z, CC_TMP1, 0, 0, 11, 0, 11, 0, "Z(nM,eI)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract424(&D, &tIA, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&D);
    /* Z(nM,eI) --> W(Mn,eI) */
    dpd_buf4_sort_axpy(&Z, CC3_HET1, qprs, 0, 11, "CC3 WMnIe (Mn,eI)", 1);
    dpd_buf4_close(&Z);

    /* Z(Nm,Ei) = <Nm|Ef> T(i,f) */
    dpd_buf4_init(&Z, CC_TMP1, 0, 0, 11, 0, 11, 0, "Z(Nm,Ei)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract424(&D, &tia, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&D);
    /* Z(Nm,Ei) --> W(mN,Ei) */
    dpd_buf4_sort_axpy(&Z, CC3_HET1, qprs, 0, 11, "CC3 WmNiE (mN,Ei)", 1);
    dpd_buf4_close(&Z);

    /* purge (mn,ei)'s before sorting */
    purge_HET1_Wmnie();

    /* also put "normal" sorted versions in CC3_HET1 */
    dpd_buf4_init(&W, CC3_HET1, 0, 2, 11, 2, 11, 0, "CC3 WMNIE (M>N,EI)");
    dpd_buf4_sort(&W, CC3_HET1, pqsr, 2, 10, "CC3 WMNIE (M>N,IE)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC3_HET1, 0, 2, 11, 2, 11, 0, "CC3 Wmnie (m>n,ei)");
    dpd_buf4_sort(&W, CC3_HET1, pqsr, 2, 10, "CC3 Wmnie (m>n,ie)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC3_HET1, 0, 0, 11, 0, 11, 0, "CC3 WMnIe (Mn,eI)");
    dpd_buf4_sort(&W, CC3_HET1, pqsr, 0, 10, "CC3 WMnIe (Mn,Ie)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC3_HET1, 0, 0, 11, 0, 11, 0, "CC3 WmNiE (mN,Ei)");
    dpd_buf4_sort(&W, CC3_HET1, pqsr, 0, 10, "CC3 WmNiE (mN,iE)");
    dpd_buf4_close(&W);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);

  }

  else if (params.ref == 2) {

    /** W(M>N,IE) <--- <MN||IE> **/
    dpd_buf4_init(&E, CC_EINTS, 0, 2, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
    dpd_buf4_sort(&E, CC3_HET1, pqsr, 2, 21, "CC3 WMNIE (M>N,EI)");
    dpd_buf4_close(&E);

    /** W(m>n,ie) <--- <mn||ie> **/
    dpd_buf4_init(&E, CC_EINTS, 0, 12, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
    dpd_buf4_sort(&E, CC3_HET1, pqsr, 12, 31, "CC3 Wmnie (m>n,ei)");
    dpd_buf4_close(&E);

    /** W(Mn,Ie) <--- <Mn|Ie> **/
    dpd_buf4_init(&E, CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    dpd_buf4_sort(&E, CC3_HET1, pqsr, 22, 25, "CC3 WMnIe (Mn,eI)");
    dpd_buf4_close(&E);

    /** W(mN,iE) <--- <mN|iE> **/
    dpd_buf4_init(&E, CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
    dpd_buf4_sort(&E, CC3_HET1, pqsr, 23, 26, "CC3 WmNiE (mN,Ei)");
    dpd_buf4_close(&E);

    /**** Term 2 ****/

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");

    /* <M>N||EF> T(I,F) --> W(M>N,EI) */
    dpd_buf4_init(&W, CC3_HET1, 0, 2, 21, 2, 21, 0, "CC3 WMNIE (M>N,EI)");
    dpd_buf4_init(&D, CC_DINTS, 0, 2, 5, 2, 5, 0, "D <IJ||AB> (I>J,AB)");
    dpd_contract424(&D, &tIA, &W, 3, 1, 0, -1, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);

    /* <m>n||ef> T(i,f) --> W(m>n,ei) */
    dpd_buf4_init(&W, CC3_HET1, 0, 12, 31, 12, 31, 0, "CC3 Wmnie (m>n,ei)");
    dpd_buf4_init(&D, CC_DINTS, 0, 12, 15, 12, 15, 0, "D <ij||ab> (i>j,ab)");
    dpd_contract424(&D, &tia, &W, 3, 1, 0, -1, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);

    /* Z(nM,eI) = <nM|eF> T(I,F) */
    dpd_buf4_init(&Z, CC_TMP1, 0, 23, 25, 23, 25, 0, "Z(nM,eI)");
    dpd_buf4_init(&D, CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    dpd_contract424(&D, &tIA, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&D);
    /* Z(nM,eI) --> W(Mn,eI) */
    dpd_buf4_sort_axpy(&Z, CC3_HET1, qprs, 22, 25, "CC3 WMnIe (Mn,eI)", 1);
    dpd_buf4_close(&Z);

    /* Z(Nm,Ei) = <Nm|Ef> T(i,f) */
    dpd_buf4_init(&Z, CC_TMP1, 0, 22, 26, 22, 26, 0, "Z(Nm,Ei)");
    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_contract424(&D, &tia, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&D);
    /* Z(Nm,Ei) --> W(mN,Ei) */
    dpd_buf4_sort_axpy(&Z, CC3_HET1, qprs, 23, 26, "CC3 WmNiE (mN,Ei)", 1);
    dpd_buf4_close(&Z);

    /* also put "normal" sorted versions in CC3_HET1 */
    dpd_buf4_init(&W, CC3_HET1, 0, 2, 21, 2, 21, 0, "CC3 WMNIE (M>N,EI)");
    dpd_buf4_sort(&W, CC3_HET1, pqsr, 2, 20, "CC3 WMNIE (M>N,IE)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC3_HET1, 0, 12, 31, 12, 31, 0, "CC3 Wmnie (m>n,ei)");
    dpd_buf4_sort(&W, CC3_HET1, pqsr, 12, 30, "CC3 Wmnie (m>n,ie)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC3_HET1, 0, 22, 25, 22, 25, 0, "CC3 WMnIe (Mn,eI)");
    dpd_buf4_sort(&W, CC3_HET1, pqsr, 22, 24, "CC3 WMnIe (Mn,Ie)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC3_HET1, 0, 23, 26, 23, 26, 0, "CC3 WmNiE (mN,Ei)");
    dpd_buf4_sort(&W, CC3_HET1, pqsr, 23, 27, "CC3 WmNiE (mN,iE)");
    dpd_buf4_close(&W);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);

  }
}
