#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void WAbEi_UHF(void)
{
  dpdfile2 Fme, T1;
  dpdbuf4 F, W, T2, B, Z, Z1, Z2, D, T, E, C;

  /** <Ei|Ab> **/
  dpd_buf4_init(&F, CC_FINTS, 0, 26, 28, 26, 28, 0, "F <Ai|Bc>");
  dpd_buf4_copy(&F, CC_HBAR, "WEiAb");
  dpd_buf4_close(&F);

  /** - F_ME t_Mi^Ab **/
  dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
  dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "FME");
  dpd_buf4_init(&W, CC_HBAR, 0, 26, 28, 26, 28, 0, "WEiAb");
  dpd_contract244(&Fme, &T2, &W, 0, 0, 0, -1.0, 1.0);
  dpd_buf4_close(&W);
  dpd_file2_close(&Fme);
  dpd_buf4_close(&T2);

  /** The next four terms are colected into Z(Ab,Ei) **/

  /** <Ab|Ef> t_i^f **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 28, 26, 28, 26, 0, "Z(Ab,Ei)");
  dpd_buf4_init(&B, CC_BINTS, 0, 28, 28, 28, 28, 0, "B <Ab|Cd>");
  dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
  dpd_contract424(&B, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&B);
  dpd_buf4_close(&Z);


  /** Z(Mb,Ei) <-- -<Mb|Ef> t_i^f **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 24, 26, 24, 26, 0, "Z(Mb,Ei)");
  dpd_buf4_init(&F, CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
  dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
  dpd_contract424(&F, &T1, &Z, 3, 1, 0, -1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&F);
  dpd_buf4_close(&Z);

  /** t_M^A Z(Mb,Ei) --> Z1(Ab,Ei) **/
  dpd_buf4_init(&Z1, CC_TMP0, 0, 28, 26, 28, 26, 0, "Z1(Ab,Ei)");
  dpd_buf4_init(&Z, CC_TMP0, 0, 24, 26, 24, 26, 0, "Z(Mb,Ei)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract244(&T1, &Z, &Z1, 0, 0, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&Z1);

  /** Z(Am,Ei) <-- -<Am|Ef> t_i^f **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 26, 26, 26, 26, 0, "Z(Am,Ei)");
  dpd_buf4_init(&F, CC_FINTS, 0, 26, 28, 26, 28, 0, "F <Ai|Bc>");
  dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
  dpd_contract424(&F, &T1, &Z, 3, 1, -1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&F);
  /** Z(Am,Ei) --> Z(mA,Ei) **/
  dpd_buf4_sort(&Z, CC_TMP0, qprs, 27, 26, "Z(mA,Ei)");
  dpd_buf4_close(&Z);

  /** t_m^b Z(mA,Ei) --> Z2(bA,Ei) **/
  dpd_buf4_init(&Z2, CC_TMP0, 0, 29, 26, 29, 26, 0, "Z(bA,Ei)");
  dpd_buf4_init(&Z, CC_TMP0, 0, 27, 26, 27, 26, 0, "Z(mA,Ei)");
  dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
  dpd_contract244(&T1, &Z, &Z2, 0, 0, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  /** Z2(bA,Ei) --> Z2(Ab,Ei) **/
  dpd_buf4_sort(&Z2, CC_TMP0, qprs, 28, 26, "Z2(Ab,Ei)");
  dpd_buf4_close(&Z2);

  /** Z1(Ab,Ei) + Z2(Ab,Ei) --> Z(Ab,Ei) **/
  dpd_buf4_init(&Z1, CC_TMP0, 0, 28, 26, 28, 26, 0, "Z1(Ab,Ei)");
  dpd_buf4_init(&Z2, CC_TMP0, 0, 28, 26, 28, 26, 0, "Z2(Ab,Ei)");
  dpd_buf4_axpy(&Z1, &Z2, 1.0);
  dpd_buf4_close(&Z1);

  dpd_buf4_init(&Z, CC_TMP0, 0, 28, 26, 28, 26, 0, "Z(Ab,Ei)");
  dpd_buf4_axpy(&Z2, &Z, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&Z2);


  dpd_buf4_init(&Z1, CC_TMP0, 0, 22, 26, 22, 26, 0, "Z(Mn,Ei)");

  /** <Mn|Ef> t_i^f --> Z1(Mn,Ei) **/
  dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
  dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
  dpd_contract424(&D, &T1, &Z1, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&D);

  /** tau_Mn^Ab Z1(Mn,Ei) --> Z(Ab,Ei) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 28, 26, 28, 26, 0, "Z(Ab,Ei)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
  dpd_contract444(&T2, &Z1, &Z, 1, 1, 1.0, 1.0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&Z);

  dpd_buf4_close(&Z1);


  /** tau_Mn^Ab <Mn|Ei> --> Z(Ab,Ei) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 28, 26, 28, 26, 0, "Z(Ab,Ei)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
  dpd_buf4_init(&E, CC_EINTS, 0, 22, 26, 22, 26, 0, "E <Ij|Ak>");
  dpd_contract444(&T2, &E, &Z, 1, 1, 1.0, 1.0);
  dpd_buf4_close(&E);
  dpd_buf4_close(&T2);

  /** Z(Ab,Ei) --> Z(Ei,Ab) --> WAbEi **/
  dpd_buf4_sort(&Z, CC_TMP0, rspq, 26, 28, "Z(Ei,Ab)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, CC_TMP0, 0, 26, 28, 26, 28, 0, "Z(Ei,Ab)");
  dpd_buf4_init(&W, CC_HBAR, 0, 26, 28, 26, 28, 0, "WEiAb");
  dpd_buf4_axpy(&Z, &W, 1.0);
  dpd_buf4_close(&W);
  dpd_buf4_close(&Z);


  /*** Remaining two terms ***/

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
