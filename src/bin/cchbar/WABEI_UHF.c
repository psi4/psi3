#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/* WABEI_UHF(): Computes all contributions to the ABEI spin case of
** the Wabei HBAR matrix elements.  The final product is stored in
** (AB,EI) ordering and is referred to on disk as "WABEI".
**
** The spin-orbital expression for the Wabei elements is:
**
** Wabei = <ab||ei> - Fme t_mi^ab + t_i^f <ab||ef>
**         - P(ab) t_m^b <am||ef>t_i^f + 1/2 tau_mn^ab <mn||ef> t_i^f
**         + 1/2 <mn||ei> tau_mn^ab - P(ab) <mb||ef> t_mi^af
**         - P(ab) t_m^a { <mb||ei> - t_ni^bf <mn||ef> }
**
** (cf. Gauss and Stanton, JCP 103, 3561-3577 (1995).)
**
** For the ABEI spin case, we evaluate these contractions with two
** target orderings, (AB,EI) and (EI,AB), depending on the term.
** After all terms have been evaluated, the (AB,EI) terms are sorted
** into (EI,AB) ordering and both groups arer added together.
**
*/

void WABEI_UHF(void)
{
  dpdfile2 Fme, T1;
  dpdbuf4 F, W, T2, B, Z, Z1, Z2, D, T, E, C;

  /** <EI||AB> **/
  dpd_buf4_init(&F, CC_FINTS, 0, 21, 7, 21, 5, 1, "F <AI|BC>");
  dpd_buf4_copy(&F, CC_HBAR, "WEIAB");
  dpd_buf4_close(&F);

  /** - F_ME t_MI^AB **/
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
  dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "FME");
  dpd_buf4_init(&W, CC_HBAR, 0, 21, 7, 21, 7, 0, "WEIAB");
  dpd_contract244(&Fme, &T2, &W, 0, 0, 0, -1.0, 1.0);
  dpd_buf4_close(&W);
  dpd_file2_close(&Fme);
  dpd_buf4_close(&T2);

  /** The next four terms are colected into Z(AB,EI) **/

  /** <AB||EF> t_I^F **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 7, 21, 7, 21, 0, "Z(AB,EI)");
  dpd_buf4_init(&B, CC_BINTS, 0, 7, 5, 5, 5, 1, "B <AB|CD>");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&B, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&B);
  dpd_buf4_close(&Z);

  /** Z(MB,EI) <-- - <MB||EF> t_I^F **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 20, 21, 20, 21, 0, "Z(MB,EI)");
  dpd_buf4_init(&F, CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&F, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&F);
  dpd_buf4_close(&Z);

  /** t_M^A Z(MB,EI) --> Z1(AB,EI) **/
  dpd_buf4_init(&Z1, CC_TMP0, 0, 5, 21, 5, 21, 0, "Z1(AB,EI)");
  dpd_buf4_init(&Z, CC_TMP0, 0, 20, 21, 20, 21, 0, "Z(MB,EI)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract244(&T1, &Z, &Z1, 0, 0, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_sort(&Z1, CC_TMP0, qprs, 5, 21, "Z2(BA,EI)");
  dpd_buf4_close(&Z1);

  /** -Z1(AB,EI) + Z2(BA,EI) --> Z(AB,EI) **/
  dpd_buf4_init(&Z1, CC_TMP0, 0, 5, 21, 5, 21, 0, "Z1(AB,EI)");
  dpd_buf4_init(&Z2, CC_TMP0, 0, 5, 21, 5, 21, 0, "Z2(BA,EI)");
  dpd_buf4_axpy(&Z1, &Z2, -1.0);
  dpd_buf4_close(&Z1);

  dpd_buf4_init(&Z, CC_TMP0, 0, 5, 21, 7, 21, 0, "Z(AB,EI)");
  dpd_buf4_axpy(&Z2, &Z, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&Z2);


  dpd_buf4_init(&Z1, CC_TMP0, 0, 2, 21, 2, 21, 0, "Z(MN,EI)");

  /** <MN||EF> t_I^F --> Z1(MN,EI) **/
  dpd_buf4_init(&D, CC_DINTS, 0, 2, 5, 2, 5, 0, "D <IJ||AB> (I>J,AB)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&D, &T1, &Z1, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&D);

  /** tau_MN^AB Z1(MN,EI) --> Z(AB,EI) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 7, 21, 7, 21, 0, "Z(AB,EI)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
  dpd_contract444(&T2, &Z1, &Z, 1, 1, 1.0, 1.0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&Z);

  dpd_buf4_close(&Z1);


  /** tau_MN^AB <MN||EI> --> Z(AB,EI) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 7, 21, 7, 21, 0, "Z(AB,EI)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
  dpd_buf4_init(&E, CC_EINTS, 0, 2, 21, 2, 21, 0, "E <IJ||KA> (I>J,AK)");
  dpd_contract444(&T2, &E, &Z, 1, 1, -1.0, 1.0);
  dpd_buf4_close(&E);
  dpd_buf4_close(&T2);

  /** Z(AB,EI) --> Z(EI,AB) --> WABEI **/
  dpd_buf4_sort(&Z, CC_TMP0, rspq, 21, 7, "Z(EI,AB)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, CC_TMP0, 0, 21, 7, 21, 7, 0, "Z(EI,AB)");
  dpd_buf4_init(&W, CC_HBAR, 0, 21, 7, 21, 7, 0, "WEIAB");
  dpd_buf4_axpy(&Z, &W, 1.0);
  dpd_buf4_close(&W);
  dpd_buf4_close(&Z);


  /** <IA||BC> --> (IB,AC) **/
  dpd_buf4_init(&F, CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
  dpd_buf4_sort(&F, CC_FINTS, prqs, 20, 5, "F <IA||BC> (IB,AC)");
  dpd_buf4_close(&F);

  /** <iA|bC> --> (ib,AC) **/
  dpd_buf4_init(&F, CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
  dpd_buf4_sort(&F, CC_FINTS, prqs, 30, 5, "F <iA|bC> (ib,AC)");
  dpd_buf4_close(&F);

  dpd_buf4_init(&Z, CC_TMP0, 0, 5, 20, 5, 20, 0, "Z(AE,IB)");

  /** <MA||FE> t_MI^FB --> Z(AE,IB) **/
  dpd_buf4_init(&F, CC_FINTS, 0, 20, 5, 20, 5, 0, "F <IA||BC> (IB,AC)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
  dpd_contract444(&F, &T2, &Z, 1, 1, 1.0, 0.0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&F);

  /** <mA|fE> t_mI^fB --> Z(AE,IB) **/
  dpd_buf4_init(&F, CC_FINTS, 0, 30, 5, 30, 5, 0, "F <iA|bC> (ib,AC)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
  dpd_contract444(&F, &T2, &Z, 1, 1, 1.0, 1.0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&F);

  dpd_buf4_close(&Z);

  dpd_buf4_init(&Z1, CC_TMP0, 0, 20, 20, 20, 20, 0, "Z1(ME,IB)");
  /** <MN||EF> t_IN^BF --> Z1(ME,IB) **/
  dpd_buf4_init(&D, CC_DINTS, 0, 20, 20, 20, 20, 0, "D <IJ||AB> (IA,JB)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
  dpd_contract444(&D, &T2, &Z1, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&D);

  /** <Mn|Ef> t_In^Bf --> Z1(ME,IB) **/
  dpd_buf4_init(&D, CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
  dpd_contract444(&D, &T2, &Z1, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&D);

  /** <IA||JB> --> (IB,JA) **/
  dpd_buf4_init(&C, CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA||JB>");
  dpd_buf4_sort(&C, CC_CINTS, psrq, 20, 20, "C <IA||JB> (IB,JA)");
  dpd_buf4_close(&C);

  /** -<MB||IE> + Z1(ME,IB) --> Z1(ME,IB) **/
  dpd_buf4_init(&C, CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA||JB> (IB,JA)");
  dpd_buf4_axpy(&C, &Z1, -1.0);
  dpd_buf4_close(&C);

  /** -t_M^A Z1(ME,IB) --> Z(AE,IB) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 5, 20, 5, 20, 0, "Z(AE,IB)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract244(&T1, &Z1, &Z, 0, 0, 0, -1.0, 1.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z1);

  /** Z(AE,IB) --> Z1(AB,EI) **/
  dpd_buf4_sort(&Z, CC_TMP0, psqr, 5, 21, "Z1(AB,EI)");
  dpd_buf4_close(&Z);

  /** Z1(AB,EI) - Z2(BA,EI) --> Z1(AB,EI) **/
  dpd_buf4_init(&Z1, CC_TMP0, 0, 5, 21, 5, 21, 0, "Z1(AB,EI)");
  dpd_buf4_sort(&Z1, CC_TMP0, qprs, 5, 21, "Z2(BA,EI)");
  dpd_buf4_init(&Z2, CC_TMP0, 0, 5, 21, 5, 21, 0, "Z2(BA,EI)");
  dpd_buf4_axpy(&Z2, &Z1, -1.0);
  dpd_buf4_close(&Z2);
  dpd_buf4_sort(&Z1, CC_TMP0, rspq, 21, 7, "Z(EI,AB)");
  dpd_buf4_close(&Z1);

  /** Z(EI,AB) --> W(EI,AB) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 21, 7, 21, 7, 0, "Z(EI,AB)");
  dpd_buf4_init(&W, CC_HBAR, 0, 21, 7, 21, 7, 0, "WEIAB");
  dpd_buf4_axpy(&Z, &W, 1.0);
  dpd_buf4_print(&W, outfile, 1);
  dpd_buf4_close(&W);
  dpd_buf4_close(&Z);
}
