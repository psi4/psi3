#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void Wabei_UHF(void)
{
  dpdfile2 Fme, T1;
  dpdbuf4 F, W, T2, B, Z, Z1, Z2, D, T, E, C;

    /** <ei||ab> **/
    dpd_buf4_init(&F, CC_FINTS, 0, 31, 17, 31, 15, 1, "F <ai|bc>");
    dpd_buf4_copy(&F, CC_HBAR, "Weiab");
    dpd_buf4_close(&F);

    /** - F_ME t_MI^AB **/
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab");
    dpd_file2_init(&Fme, CC_OEI, 0, 2, 3, "Fme");
    dpd_buf4_init(&W, CC_HBAR, 0, 31, 17, 31, 17, 0, "Weiab");
    dpd_contract244(&Fme, &T2, &W, 0, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_file2_close(&Fme);
    dpd_buf4_close(&T2);

    /** The next four terms are colected into Z(AB,EI) **/

    /** <AB||EF> t_I^F **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 17, 31, 17, 31, 0, "Z(ab,ei)");
    dpd_buf4_init(&B, CC_BINTS, 0, 17, 15, 15, 15, 1, "B <ab|cd>");
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_contract424(&B, &T1, &Z, 3, 1, 0, 1.0, 0.0);
    dpd_file2_close(&T1);
    dpd_buf4_close(&B);
    dpd_buf4_close(&Z);

    /** Z(MB,EI) <-- - <MB||EF> t_I^F **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 30, 31, 30, 31, 0, "Z(mb,ei)");
    dpd_buf4_init(&F, CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_contract424(&F, &T1, &Z, 3, 1, 0, 1.0, 0.0);
    dpd_file2_close(&T1);
    dpd_buf4_close(&F);
    dpd_buf4_close(&Z);

    /** t_M^A Z(MB,EI) --> Z1(AB,EI) **/
    dpd_buf4_init(&Z1, CC_TMP0, 0, 15, 31, 15, 31, 0, "Z1(ab,ei)");
    dpd_buf4_init(&Z, CC_TMP0, 0, 30, 31, 30, 31, 0, "Z(mb,ei)");
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_contract244(&T1, &Z, &Z1, 0, 0, 0, 1.0, 0.0);
    dpd_file2_close(&T1);
    dpd_buf4_close(&Z);
    dpd_buf4_sort(&Z1, CC_TMP0, qprs, 15, 31, "Z2(ba,ei)");
    dpd_buf4_close(&Z1);

    /** -Z1(AB,EI) + Z2(BA,EI) --> Z(AB,EI) **/
    dpd_buf4_init(&Z1, CC_TMP0, 0, 15, 31, 15, 31, 0, "Z1(ab,ei)");
    dpd_buf4_init(&Z2, CC_TMP0, 0, 15, 31, 15, 31, 0, "Z2(ba,ei)");
    dpd_buf4_axpy(&Z1, &Z2, -1.0);
    dpd_buf4_close(&Z1);

    dpd_buf4_init(&Z, CC_TMP0, 0, 15, 31, 17, 31, 0, "Z(ab,ei)");
    dpd_buf4_axpy(&Z2, &Z, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&Z2);


    dpd_buf4_init(&Z1, CC_TMP0, 0, 12, 31, 12, 31, 0, "Z(mn,ei)");

    /** <MN||EF> t_I^F --> Z1(MN,EI) **/
    dpd_buf4_init(&D, CC_DINTS, 0, 12, 15, 12, 15, 0, "D <ij||ab> (i>j,ab)");
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_contract424(&D, &T1, &Z1, 3, 1, 0, 1.0, 0.0);
    dpd_file2_close(&T1);
    dpd_buf4_close(&D);

    /** tau_MN^AB Z1(MN,EI) --> Z(AB,EI) **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 17, 31, 17, 31, 0, "Z(ab,ei)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
    dpd_contract444(&T2, &Z1, &Z, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&Z);

    dpd_buf4_close(&Z1);


    /** tau_MN^AB <MN||EI> --> Z(AB,EI) **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 17, 31, 17, 31, 0, "Z(ab,ei)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
    dpd_buf4_init(&E, CC_EINTS, 0, 12, 31, 12, 31, 0, "E <ij||ka> (i>j,ak)");
    dpd_contract444(&T2, &E, &Z, 1, 1, -1.0, 1.0);
    dpd_buf4_close(&E);
    dpd_buf4_close(&T2);

    /** Z(AB,EI) --> Z(EI,AB) --> WABEI **/
    dpd_buf4_sort(&Z, CC_TMP0, rspq, 31, 17, "Z(ei,ab)");
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, CC_TMP0, 0, 31, 17, 31, 17, 0, "Z(ei,ab)");
    dpd_buf4_init(&W, CC_HBAR, 0, 31, 17, 31, 17, 0, "Weiab");
    dpd_buf4_axpy(&Z, &W, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Z);


    /** <IA||BC> --> (IB,AC) **/
    dpd_buf4_init(&F, CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
    dpd_buf4_sort(&F, CC_FINTS, prqs, 30, 15, "F <ia||bc> (ib,ac)");
    dpd_buf4_close(&F);

    /** <iA|bC> --> (ib,AC) **/
    dpd_buf4_init(&F, CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    dpd_buf4_sort(&F, CC_FINTS, prqs, 20, 15, "F <Ia|Bc> (IB,ac)");
    dpd_buf4_close(&F);

    dpd_buf4_init(&Z, CC_TMP0, 0, 15, 30, 15, 30, 0, "Z(ae,ib)");

    /** <MA||FE> t_MI^FB --> Z(AE,IB) **/
    dpd_buf4_init(&F, CC_FINTS, 0, 30, 15, 30, 15, 0, "F <ia||bc> (ib,ac)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
    dpd_contract444(&F, &T2, &Z, 1, 1, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&F);

    /** <mA|fE> t_mI^fB --> Z(AE,IB) **/
    dpd_buf4_init(&F, CC_FINTS, 0, 20, 15, 20, 15, 0, "F <Ia|Bc> (IB,ac)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    dpd_contract444(&F, &T2, &Z, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&F);

    dpd_buf4_close(&Z);

    dpd_buf4_init(&Z1, CC_TMP0, 0, 30, 30, 30, 30, 0, "Z1(me,ib)");
    /** <MN||EF> t_IN^BF --> Z1(ME,IB) **/
    dpd_buf4_init(&D, CC_DINTS, 0, 30, 30, 30, 30, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
    dpd_contract444(&D, &T2, &Z1, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&D);

    /** <Mn|Ef> t_In^Bf --> Z1(ME,IB) **/
    dpd_buf4_init(&D, CC_DINTS, 0, 30, 20, 30, 20, 0, "D <Ij|Ab> (ia,JB)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
    dpd_contract444(&D, &T2, &Z1, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&D);

    /** <IA||JB> --> (IB,JA) **/
    dpd_buf4_init(&C, CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia||jb>");
    dpd_buf4_sort(&C, CC_CINTS, psrq, 30, 30, "C <ia||jb> (ib,ja)");
    dpd_buf4_close(&C);

    /** -<MB||IE> + Z1(ME,IB) --> Z1(ME,IB) **/
    dpd_buf4_init(&C, CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia||jb> (ib,ja)");
    dpd_buf4_axpy(&C, &Z1, -1.0);
    dpd_buf4_close(&C);

    /** -t_M^A Z1(ME,IB) --> Z(AE,IB) **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 15, 30, 15, 30, 0, "Z(ae,ib)");
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_contract244(&T1, &Z1, &Z, 0, 0, 0, -1.0, 1.0);
    dpd_file2_close(&T1);
    dpd_buf4_close(&Z1);

    /** Z(AE,IB) --> Z1(AB,EI) **/
    dpd_buf4_sort(&Z, CC_TMP0, psqr, 15, 31, "Z1(ab,ei)");
    dpd_buf4_close(&Z);

    /** Z1(AB,EI) - Z2(BA,EI) --> Z1(AB,EI) **/
    dpd_buf4_init(&Z1, CC_TMP0, 0, 15, 31, 15, 31, 0, "Z1(ab,ei)");
    dpd_buf4_sort(&Z1, CC_TMP0, qprs, 15, 31, "Z2(ba,ei)");
    dpd_buf4_init(&Z2, CC_TMP0, 0, 15, 31, 15, 31, 0, "Z2(ba,ei)");
    dpd_buf4_axpy(&Z2, &Z1, -1.0);
    dpd_buf4_close(&Z2);
    dpd_buf4_sort(&Z1, CC_TMP0, rspq, 31, 17, "Z(ei,ab)");
    dpd_buf4_close(&Z1);

    /** Z(EI,AB) --> W(EI,AB) **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 31, 17, 31, 17, 0, "Z(ei,ab)");
    dpd_buf4_init(&W, CC_HBAR, 0, 31, 17, 31, 17, 0, "Weiab");
    dpd_buf4_axpy(&Z, &W, 1.0);
    dpd_buf4_print(&W, outfile, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Z);
}
