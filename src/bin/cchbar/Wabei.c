#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void Wabei_build(void)
{
  struct oe_dpdfile Fme, T1;
  struct dpdbuf F, W, T2, B, Z, Z1, Z2, D, T, E, C;
  
  dpd_buf_init(&F, CC_FINTS, 10, 7, 10, 5, 1, "F <ia|bc>", 0, outfile);
  /** <EI||AB> **/
  dpd_swap12(&F, CC_HBAR, 11, 7, "WEIAB", 0, outfile);
  /** <ei||ab> **/
  dpd_swap12(&F, CC_HBAR, 11, 7, "Weiab", 0, outfile);
  dpd_buf_close(&F);

  dpd_buf_init(&W, CC_HBAR, 11, 7, 11, 7, 0, "WEIAB", 0, outfile);
  dpd_scm(&W, -1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_init(&W, CC_HBAR, 11, 7, 11, 7, 0, "Weiab", 0, outfile);
  dpd_scm(&W, -1.0, 0, outfile);
  dpd_buf_close(&W);
  
  dpd_buf_init(&F, CC_FINTS, 10, 5, 10, 5, 0, "F <ia|bc>", 0, outfile);
  /** <iE|bA> **/
  dpd_swap12(&F, CC_TMP0, 11, 5, "W(Ei,bA)", 0, outfile);
  dpd_buf_close(&F);
  dpd_buf_init(&W, CC_TMP0, 11, 5, 11, 5, 0, "W(Ei,bA)", 0, outfile);
  dpd_swap34(&W, CC_HBAR, 11, 5, "WEiAb", 0, outfile);
  dpd_buf_close(&W);
  /** <Ie|Ba> **/
  dpd_buf_init(&W, CC_HBAR, 11, 5, 11, 5, 0, "WEiAb", 0, outfile);
  dpd_copy(&W, CC_HBAR, "WeIaB", 0, outfile);
  dpd_buf_close(&W);

  /** - F_ME t_MI^AB **/
  dpd_buf_init(&T2, CC_TAMPS, 0, 7, 2, 7, 0, "tIJAB", 0, outfile);
  dpd_oe_file_init(&Fme, CC_OEI, 0, 1, "FME", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 11, 7, 11, 7, 0, "WEIAB", 0, outfile);
  dpd_contract212(&Fme, &T2, &W, 0, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_oe_file_close(&Fme);
  dpd_buf_close(&T2);

  /** - F_me t_mi^ab **/
  dpd_buf_init(&T2, CC_TAMPS, 0, 7, 2, 7, 0, "tijab", 0, outfile);
  dpd_oe_file_init(&Fme, CC_OEI, 0, 1, "Fme", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 11, 7, 11, 7, 0, "Weiab", 0, outfile);
  dpd_contract212(&Fme, &T2, &W, 0, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_oe_file_close(&Fme);
  dpd_buf_close(&T2);

  /** - F_ME t_Mi^Ab **/
  dpd_buf_init(&T2, CC_TAMPS, 0, 5, 0, 5, 0, "tIjAb", 0, outfile);
  dpd_oe_file_init(&Fme, CC_OEI, 0, 1, "FME", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 11, 5, 11, 5, 0, "WEiAb", 0, outfile);
  dpd_contract212(&Fme, &T2, &W, 0, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_oe_file_close(&Fme);
  dpd_buf_close(&T2);

  /** - F_me t_mI^aB **/
  dpd_buf_init(&T2, CC_TAMPS, 0, 5, 0, 5, 0, "tiJaB", 0, outfile);
  dpd_oe_file_init(&Fme, CC_OEI, 0, 1, "Fme", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 11, 5, 11, 5, 0, "WeIaB", 0, outfile);
  dpd_contract212(&Fme, &T2, &W, 0, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_oe_file_close(&Fme);
  dpd_buf_close(&T2);

  /** <AB||EF> t_I^F **/
  dpd_buf_init(&Z, CC_TMP0, 7, 11, 7, 11, 0, "Z(AB,EI)", 0, outfile);
  dpd_buf_init(&B, CC_BINTS, 7, 5, 5, 5, 1, "B <ab|cd>", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_contract221(&B, &T1, &Z, 3, 1, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&B);
  dpd_buf_sort(&Z, CC_TMP1, rspq, 11, 7, "Z(EI,AB)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&Z, CC_TMP1, 11, 7, 11, 7, 0, "Z(EI,AB)", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 11, 7, 11, 7, 0, "WEIAB", 0, outfile);
  dpd_axpy(&Z, &W, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&Z);

  /** <ab||ef> t_i^f **/
  dpd_buf_init(&Z, CC_TMP0, 7, 11, 7, 11, 0, "Z(ab,ei)", 0, outfile);
  dpd_buf_init(&B, CC_BINTS, 7, 5, 5, 5, 1, "B <ab|cd>", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_contract221(&B, &T1, &Z, 3, 1, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&B);
  dpd_buf_sort(&Z, CC_TMP1, rspq, 11, 7, "Z(ei,ab)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&Z, CC_TMP1, 11, 7, 11, 7, 0, "Z(ei,ab)", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 11, 7, 11, 7, 0, "Weiab", 0, outfile);
  dpd_axpy(&Z, &W, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&Z);

  /** - <fE|bA> t_i^f **/
  dpd_buf_init(&Z, CC_TMP0, 5, 11, 5, 11, 0, "Z(Ab,Ei)", 0, outfile);
  dpd_buf_init(&B, CC_BINTS, 5, 5, 5, 5, 0, "B <ab|cd>", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_contract221(&B, &T1, &Z, 3, 1, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&B);
  dpd_buf_sort(&Z, CC_TMP1, rspq, 11, 5, "Z(Ei,Ab)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&Z, CC_TMP1, 11, 5, 11, 5, 0, "Z(Ei,Ab)", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 11, 5, 11, 5, 0, "WEiAb", 0, outfile);
  dpd_axpy(&Z, &W, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&Z);

  /** - <Fe|Ba> t_I^F **/
  dpd_buf_init(&Z, CC_TMP0, 5, 11, 5, 11, 0, "Z(aB,eI)", 0, outfile);
  dpd_buf_init(&B, CC_BINTS, 5, 5, 5, 5, 0, "B <ab|cd>", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_contract221(&B, &T1, &Z, 3, 1, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&B);
  dpd_buf_sort(&Z, CC_TMP1, rspq, 11, 5, "Z(eI,aB)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&Z, CC_TMP1, 11, 5, 11, 5, 0, "Z(eI,aB)", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 11, 5, 11, 5, 0, "WeIaB", 0, outfile);
  dpd_axpy(&Z, &W, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&Z);

  /** Prepare intermediates for second Wabef contribution to Wabei **/
  /** Z(MA,EI) <-- <MA||EF> t_I^F **/
  dpd_buf_init(&Z, CC_TMP0, 10, 11, 10, 11, 0, "Z(MA,EI)", 0, outfile);
  dpd_buf_init(&F, CC_FINTS, 10, 5, 10, 5, 1, "F <ia|bc>", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_contract221(&F, &T1, &Z, 3, 1, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&F);
  dpd_buf_close(&Z);

  /** t_M^B Z(MA,EI) --> Z'(BA,EI) --> Z1(AB,EI) **/
  dpd_buf_init(&Z1, CC_TMP1, 5, 11, 5, 11, 0, "Z(BA,EI)", 0, outfile);
  dpd_buf_init(&Z, CC_TMP0, 10, 11, 10, 11, 0, "Z(MA,EI)", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_contract212(&T1, &Z, &Z1, 0, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&Z);
  dpd_swap12(&Z1, CC_TMP2, 5, 11, "Z(AB,EI)", 0, outfile);
  dpd_buf_close(&Z1);

  /** t_M^A Z(MB,EI) --> Z2(AB,EI) **/
  dpd_buf_init(&Z2, CC_TMP1, 5, 11, 5, 11, 0, "Z(AB,EI)", 0, outfile);
  dpd_buf_init(&Z, CC_TMP0, 10, 11, 10, 11, 0, "Z(MA,EI)", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_contract212(&T1, &Z, &Z2, 0, 0, 0, -1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&Z);

  /** Z1(AB,EI) + Z2(AB,EI) --> W(AB,EI) --> W(EI,AB) **/
  dpd_buf_init(&Z1, CC_TMP2, 5, 11, 5, 11, 0, "Z(AB,EI)", 0, outfile);
  dpd_axpy(&Z1, &Z2, 1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_swapbk(&Z2, CC_TMP0, 11, 5, "Z(EI,AB)", 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_init(&Z2, CC_TMP0, 11, 5, 11, 5, 0, "Z(EI,AB)", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 11, 5, 11, 7, 0, "WEIAB", 0, outfile);
  dpd_axpy(&Z2, &W, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&Z2);

  
  /** Z(ma,ei) <-- <ma||ef> t_i^f **/
  dpd_buf_init(&Z, CC_TMP0, 10, 11, 10, 11, 0, "Z(ma,ei)", 0, outfile);
  dpd_buf_init(&F, CC_FINTS, 10, 5, 10, 5, 1, "F <ia|bc>", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_contract221(&F, &T1, &Z, 3, 1, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&F);
  dpd_buf_close(&Z);

  /** t_m^b Z(ma,ei) --> Z'(ba,ei) --> Z1(ab,ei) **/
  dpd_buf_init(&Z1, CC_TMP1, 5, 11, 5, 11, 0, "Z(ba,ei)", 0, outfile);
  dpd_buf_init(&Z, CC_TMP0, 10, 11, 10, 11, 0, "Z(ma,ei)", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_contract212(&T1, &Z, &Z1, 0, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&Z);
  dpd_swap12(&Z1, CC_TMP2, 5, 11, "Z(ab,ei)", 0, outfile);
  dpd_buf_close(&Z1);

  /** t_m^a Z(mb,ei) --> Z2(ab,ei) **/
  dpd_buf_init(&Z2, CC_TMP1, 5, 11, 5, 11, 0, "Z(ab,ei)", 0, outfile);
  dpd_buf_init(&Z, CC_TMP0, 10, 11, 10, 11, 0, "Z(ma,ei)", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_contract212(&T1, &Z, &Z2, 0, 0, 0, -1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&Z);

  /** Z1(ab,ei) + Z2(ab,ei) --> W(ab,ei) --> W(ie,ab) **/
  dpd_buf_init(&Z1, CC_TMP2, 5, 11, 5, 11, 0, "Z(ab,ei)", 0, outfile);
  dpd_axpy(&Z1, &Z2, 1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_swapbk(&Z2, CC_TMP0, 11, 5, "Z(ei,ab)", 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_init(&Z2, CC_TMP0, 11, 5, 11, 5, 0, "Z(ei,ab)", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 11, 5, 11, 7, 0, "Weiab", 0, outfile);
  dpd_axpy(&Z2, &W, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&Z2);

  /** t_i^f <mA|fE> --> Z(mA,Ei) **/
  dpd_buf_init(&Z, CC_TMP0, 10, 10, 10, 10, 0, "Z(mA,iE)", 0, outfile);
  dpd_buf_init(&F, CC_FINTS, 10, 5, 10, 5, 0, "F <ia|bc>", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_contract212(&T1, &F, &Z, 1, 2, 1, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&F);
  dpd_swap34(&Z, CC_TMP1, 10, 11, "Z(mA,Ei)", 0, outfile);
  dpd_buf_close(&Z);

  /** <Mb|Ef> t_i^f --> Z(Mb,Ei) **/
  dpd_buf_init(&Z, CC_TMP0, 10, 11, 10, 11, 0, "Z(Mb,Ei)", 0, outfile);
  dpd_buf_init(&F, CC_FINTS, 10, 5, 10, 5, 0, "F <ia|bc>", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_contract221(&F, &T1, &Z, 3, 1, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&F);
  dpd_buf_close(&Z);

  /** - T_M^A Z(Mb,Ei) --> Z(Ab,Ei) **/
  dpd_buf_init(&Z1, CC_TMP2, 5, 11, 5, 11, 0, "Z(Ab,Ei)", 0, outfile);
  dpd_buf_init(&Z, CC_TMP0, 10, 11, 10, 11, 0, "Z(Mb,Ei)", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_contract212(&T1, &Z, &Z1, 0, 0, 0, -1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&Z);
  dpd_buf_close(&Z1);

  /** - t_m^b Z(mA,Ei) --> Z1(Ab,Ei) **/
  dpd_buf_init(&Z1, CC_TMP0, 5, 11, 5, 11, 0, "Z1(bA,Ei)", 0, outfile);
  dpd_buf_init(&Z, CC_TMP1, 10, 11, 10, 11, 0, "Z(mA,Ei)", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_contract212(&T1, &Z, &Z1, 0, 0, 0, -1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&Z);
  dpd_swap12(&Z1, CC_TMP1, 5, 11, "Z(Ab,Ei)", 0, outfile);
  dpd_buf_close(&Z1);

  /** Z1(Ab,Ei) + Z2(Ab,Ei) --> W(Ab,Ei) **/
  dpd_buf_init(&Z1, CC_TMP1, 5, 11, 5, 11, 0, "Z(Ab,Ei)", 0, outfile);
  dpd_buf_init(&Z2, CC_TMP2, 5, 11, 5, 11, 0, "Z(Ab,Ei)", 0, outfile);
  dpd_axpy(&Z1, &Z2, 1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_swapbk(&Z2, CC_TMP0, 11, 5, "Z(Ei,Ab)", 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_init(&Z, CC_TMP0, 11, 5, 11, 5, 0, "Z(Ei,Ab)", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 11, 5, 11, 5, 0, "WEiAb", 0, outfile);
  dpd_axpy(&Z, &W, 1.0, 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&W);

  /** t_I^F <Ma|Fe> --> Z(Ma,eI) **/
  dpd_buf_init(&Z, CC_TMP0, 10, 10, 10, 10, 0, "Z(Ma,Ie)", 0, outfile);
  dpd_buf_init(&F, CC_FINTS, 10, 5, 10, 5, 0, "F <ia|bc>", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_contract212(&T1, &F, &Z, 1, 2, 1, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&F);
  dpd_swap34(&Z, CC_TMP1, 10, 11, "Z(Ma,eI)", 0, outfile);
  dpd_buf_close(&Z);

  /** <mB|eF> t_I^F --> Z(mB,eI) **/
  dpd_buf_init(&Z, CC_TMP0, 10, 11, 10, 11, 0, "Z(mB,eI)", 0, outfile);
  dpd_buf_init(&F, CC_FINTS, 10, 5, 10, 5, 0, "F <ia|bc>", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_contract221(&F, &T1, &Z, 3, 1, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&F);
  dpd_buf_close(&Z);

  /** t_m^a Z(mB,eI) --> Z(aB,eI) **/
  dpd_buf_init(&Z1, CC_TMP2, 5, 11, 5, 11, 0, "Z(aB,eI)", 0, outfile);
  dpd_buf_init(&Z, CC_TMP0, 10, 11, 10, 11, 0, "Z(mB,eI)", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_contract212(&T1, &Z, &Z1, 0, 0, 0, -1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&Z);
  dpd_buf_close(&Z1);

  /** t_M^B Z(Ma,eI) --> Z1(aB,eI) **/
  dpd_buf_init(&Z1, CC_TMP0, 5, 11, 5, 11, 0, "Z1(Ba,eI)", 0, outfile);
  dpd_buf_init(&Z, CC_TMP1, 10, 11, 10, 11, 0, "Z(Ma,eI)", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_contract212(&T1, &Z, &Z1, 0, 0, 0, -1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&Z);
  dpd_swap12(&Z1, CC_TMP1, 5, 11, "Z(aB,eI)", 0, outfile);
  dpd_buf_close(&Z1);

  /** Z1(aB,eI) + Z2(aB,eI) --> W(aB,eI) **/
  dpd_buf_init(&Z1, CC_TMP1, 5, 11, 5, 11, 0, "Z(aB,eI)", 0, outfile);
  dpd_buf_init(&Z2, CC_TMP2, 5, 11, 5, 11, 0, "Z(aB,eI)", 0, outfile);
  dpd_axpy(&Z1, &Z2, 1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_swapbk(&Z2, CC_TMP0, 11, 5, "Z(eI,aB)", 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_init(&Z, CC_TMP0, 11, 5, 11, 5, 0, "Z(eI,aB)", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 11, 5, 11, 5, 0, "WeIaB", 0, outfile);
  dpd_axpy(&Z, &W, 1.0, 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&W);

  /** Final term of Wabef contribution to Wabei **/

  /** t_I^F <MN||EF> --> Z(MN,EI) **/
  dpd_buf_init(&Z, CC_TMP0, 2, 11, 2, 11, 0, "Z(MN,EI)", 0, outfile);
  dpd_buf_init(&D, CC_DINTS, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)",
	       0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_contract221(&D, &T1, &Z, 3, 1, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&D);
  /** Z(MN,EI) Tau(MN,AB) --> W(EI,AB) **/
  dpd_buf_init(&T, CC_TAMPS, 2, 7, 2, 7, 0, "tauIJAB", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 11, 7, 11, 7, 0, "WEIAB", 0, outfile);
  dpd_contract222(&Z, &T, &W, 1, 1, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&T);
  dpd_buf_close(&Z);

  /** t_i^f <mn||ef> --> Z(mn,ei) **/
  dpd_buf_init(&Z, CC_TMP0, 2, 11, 2, 11, 0, "Z(mn,ei)", 0, outfile);
  dpd_buf_init(&D, CC_DINTS, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)",
	       0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_contract221(&D, &T1, &Z, 3, 1, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&D);
  /** Z(mn,ei) Tau(mn,ab) --> W(ei,ab) **/
  dpd_buf_init(&T, CC_TAMPS, 2, 7, 2, 7, 0, "tauijab", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 11, 7, 11, 7, 0, "Weiab", 0, outfile);
  dpd_contract222(&Z, &T, &W, 1, 1, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&T);
  dpd_buf_close(&Z);

  /** t_i^f <Mn|Ef> --> Z(Mn,Ei) **/
  dpd_buf_init(&Z, CC_TMP0, 0, 11, 0, 11, 0, "Z(Mn,Ei)", 0, outfile);
  dpd_buf_init(&D, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_contract221(&D, &T1, &Z, 3, 1, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&D);
  /** Z(Mn,Ei) Tau(Mn,Ab) --> W(Ei,Ab) **/
  dpd_buf_init(&T, CC_TAMPS, 0, 5, 0, 5, 0, "tauIjAb", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 11, 5, 11, 5, 0, "WEiAb", 0, outfile);
  dpd_contract222(&Z, &T, &W, 1, 1, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&T);
  dpd_buf_close(&Z);

  /** t_I^F <mN|eF> --> Z(mN,eI) **/
  dpd_buf_init(&Z, CC_TMP0, 0, 11, 0, 11, 0, "Z(mN,eI)", 0, outfile);
  dpd_buf_init(&D, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_contract221(&D, &T1, &Z, 3, 1, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&D);
  /** Z(mN,eI) Tau(mN,aB) --> W(eI,aB) **/
  dpd_buf_init(&T, CC_TAMPS, 0, 5, 0, 5, 0, "tauiJaB", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 11, 5, 11, 5, 0, "WeIaB", 0, outfile);
  dpd_contract222(&Z, &T, &W, 1, 1, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&T);
  dpd_buf_close(&Z);


  /** <MN||EI> Tau(MN,AB) --> W(EI,AB) **/
  dpd_buf_init(&E, CC_EINTS, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)",
	       0, outfile);
  dpd_buf_init(&T, CC_TAMPS, 2, 7, 2, 7, 0, "tauIJAB", 0, outfile);
  dpd_buf_init(&Z, CC_TMP0, 10, 7, 10, 7, 0, "Z(IE,AB)", 0, outfile);
  dpd_contract222(&E, &T, &Z, 1, 1, -1.0, 0.0, 0, outfile);
  dpd_swap12(&Z, CC_TMP1, 11, 7, "Z(EI,AB)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&T);
  dpd_buf_close(&E);
  dpd_buf_init(&Z, CC_TMP1, 11, 7, 11, 7, 0, "Z(EI,AB)", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 11, 7, 11, 7, 0, "WEIAB", 0, outfile);
  dpd_axpy(&Z, &W, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&Z);

  /** <mn||ei> Tau(mn,ab) --> W(ei,ab) **/
  dpd_buf_init(&E, CC_EINTS, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)",
	       0, outfile);
  dpd_buf_init(&T, CC_TAMPS, 2, 7, 2, 7, 0, "tauijab", 0, outfile);
  dpd_buf_init(&Z, CC_TMP0, 10, 7, 10, 7, 0, "Z(ie,ab)", 0, outfile);
  dpd_contract222(&E, &T, &Z, 1, 1, -1.0, 0.0, 0, outfile);
  dpd_swap12(&Z, CC_TMP1, 11, 7, "Z(ei,ab)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&T);
  dpd_buf_close(&E);
  dpd_buf_init(&Z, CC_TMP1, 11, 7, 11, 7, 0, "Z(ei,ab)", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 11, 7, 11, 7, 0, "Weiab", 0, outfile);
  dpd_axpy(&Z, &W, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&Z);

  /** <Mn|Ei> Tau(Mn,Ab) --> W(Ei,Ab) **/
  dpd_buf_init(&E, CC_EINTS, 0, 10, 0, 10, 0, "E <ij|ka>", 0, outfile);
  dpd_buf_init(&T, CC_TAMPS, 0, 5, 0, 5, 0, "tauiJaB", 0, outfile);
  dpd_buf_init(&Z, CC_TMP0, 10, 5, 10, 5, 0, "Z(iE,bA)", 0, outfile);
  dpd_contract222(&E, &T, &Z, 1, 1, 1.0, 0.0, 0, outfile);
  dpd_swap12(&Z, CC_TMP1, 11, 5, "Z(Ei,bA)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&T);
  dpd_buf_close(&E);
  dpd_buf_init(&Z, CC_TMP1, 11, 5, 11, 5, 0, "Z(Ei,bA)", 0, outfile);
  dpd_swap34(&Z, CC_TMP0, 11, 5, "Z(Ei,Ab)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&Z, CC_TMP0, 11, 5, 11, 5, 0, "Z(Ei,Ab)", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 11, 5, 11, 5, 0, "WEiAb", 0, outfile);
  dpd_axpy(&Z, &W, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&Z);

  /** <mN|eI> Tau(mN,aB) --> W(eI,aB) **/
  dpd_buf_init(&E, CC_EINTS, 0, 10, 0, 10, 0, "E <ij|ka>", 0, outfile);
  dpd_buf_init(&T, CC_TAMPS, 0, 5, 0, 5, 0, "tauIjAb", 0, outfile);
  dpd_buf_init(&Z, CC_TMP0, 10, 5, 10, 5, 0, "Z(Ie,Ba)", 0, outfile);
  dpd_contract222(&E, &T, &Z, 1, 1, 1.0, 0.0, 0, outfile);
  dpd_swap12(&Z, CC_TMP1, 11, 5, "Z(eI,Ba)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&T);
  dpd_buf_close(&E);
  dpd_buf_init(&Z, CC_TMP1, 11, 5, 11, 5, 0, "Z(eI,Ba)", 0, outfile);
  dpd_swap34(&Z, CC_TMP0, 11, 5, "Z(eI,aB)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&Z, CC_TMP0, 11, 5, 11, 5, 0, "Z(eI,aB)", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 11, 5, 11, 5, 0, "WeIaB", 0, outfile);
  dpd_axpy(&Z, &W, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&Z);

  /** <MB||EF> t_IM^AF + <MA||FE> t_IM^BF --> W(EI,AB) **/
  dpd_buf_init(&F, CC_FINTS, 10, 5, 10, 5, 1, "F <ia|bc>", 0, outfile);
  dpd_swap23(&F, CC_TMP0, 10, 5, "F(MF,AE)", 0, outfile);
  dpd_buf_close(&F);
  dpd_buf_init(&F, CC_TMP0, 10, 5, 10, 5, 0, "F(MF,AE)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tIAJB", 0, outfile);
  dpd_buf_init(&Z, CC_TMP1, 10, 5, 10, 5, 0, "Z(IB,AE)", 0, outfile);
  dpd_contract222(&T2, &F, &Z, 0, 1, 1.0, 0.0, 0, outfile);
  dpd_swap24(&Z, CC_TMP0, 10, 5, "Z(IE,AB)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&T2);
  dpd_buf_close(&F);
  dpd_buf_init(&Z, CC_TMP0, 10, 5, 10, 5, 0, "Z(IE,AB)", 0, outfile);
  dpd_swap12(&Z, CC_TMP1, 11, 5, "Z(EI,AB)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&Z1, CC_TMP1, 11, 5, 11, 5, 0, "Z(EI,AB)", 0, outfile);
  dpd_swap34(&Z1, CC_TMP0, 11, 5, "Z(EI,BA)", 0, outfile);
  dpd_buf_init(&Z2, CC_TMP0, 11, 5, 11, 5, 0, "Z(EI,BA)", 0, outfile);
  dpd_axpy(&Z2, &Z1, -1.0, 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_init(&W, CC_HBAR, 11, 5, 11, 7, 0, "WEIAB", 0, outfile);
  dpd_axpy(&Z1, &W, 1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_close(&W);

  /** -<mB|fE> t_Im^Af + <mA|fE> t_Im^Bf --> W(EI,AB) **/
  dpd_buf_init(&F, CC_FINTS, 10, 5, 10, 5, 0, "F <ia|bc>", 0, outfile);
  dpd_swap23(&F, CC_TMP0, 10, 5, "F(mf,AE)", 0, outfile);
  dpd_buf_close(&F);
  dpd_buf_init(&F, CC_TMP0, 10, 5, 10, 5, 0, "F(mf,AE)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tIAjb", 0, outfile);
  dpd_buf_init(&Z, CC_TMP1, 10, 5, 10, 5, 0, "Z(IB,AE)", 0, outfile);
  dpd_contract222(&T2, &F, &Z, 0, 1, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&F);
  dpd_swap24(&Z, CC_TMP0, 10, 5, "Z(IE,AB)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&Z, CC_TMP0, 10, 5, 10, 5, 0, "Z(IE,AB)", 0, outfile);
  dpd_swap12(&Z, CC_TMP1, 11, 5, "Z(EI,AB)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&Z1, CC_TMP1, 11, 5, 11, 5, 0, "Z(EI,AB)", 0, outfile);
  dpd_swap34(&Z1, CC_TMP0, 11, 5, "Z(EI,BA)", 0, outfile);
  dpd_buf_init(&Z2, CC_TMP0, 11, 5, 11, 5, 0, "Z(EI,BA)", 0, outfile);
  dpd_axpy(&Z2, &Z1, -1.0, 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_init(&W, CC_HBAR, 11, 5, 11, 7, 0, "WEIAB", 0, outfile);
  dpd_axpy(&Z1, &W, 1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_close(&W);

  /** <mb||ef> t_im^af + <ma||fe> t_im^bf --> W(ei,ab) **/
  dpd_buf_init(&F, CC_FINTS, 10, 5, 10, 5, 1, "F <ia|bc>", 0, outfile);
  dpd_swap23(&F, CC_TMP0, 10, 5, "F(mf,ae)", 0, outfile);
  dpd_buf_close(&F);
  dpd_buf_init(&F, CC_TMP0, 10, 5, 10, 5, 0, "F(mf,ae)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tiajb", 0, outfile);
  dpd_buf_init(&Z, CC_TMP1, 10, 5, 10, 5, 0, "Z(ib,ae)", 0, outfile);
  dpd_contract222(&T2, &F, &Z, 0, 1, 1.0, 0.0, 0, outfile);
  dpd_swap24(&Z, CC_TMP0, 10, 5, "Z(ie,ab)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&T2);
  dpd_buf_close(&F);
  dpd_buf_init(&Z, CC_TMP0, 10, 5, 10, 5, 0, "Z(ie,ab)", 0, outfile);
  dpd_swap12(&Z, CC_TMP1, 11, 5, "Z(ei,ab)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&Z1, CC_TMP1, 11, 5, 11, 5, 0, "Z(ei,ab)", 0, outfile);
  dpd_swap34(&Z1, CC_TMP0, 11, 5, "Z(ei,ba)", 0, outfile);
  dpd_buf_init(&Z2, CC_TMP0, 11, 5, 11, 5, 0, "Z(ei,ba)", 0, outfile);
  dpd_axpy(&Z2, &Z1, -1.0, 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_init(&W, CC_HBAR, 11, 5, 11, 7, 0, "Weiab", 0, outfile);
  dpd_axpy(&Z1, &W, 1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_close(&W);

  /** -<Mb|Fe> t_iM^aF + <Ma|Fe> t_iM^bF --> W(ei,ab) **/
  dpd_buf_init(&F, CC_FINTS, 10, 5, 10, 5, 0, "F <ia|bc>", 0, outfile);
  dpd_swap23(&F, CC_TMP0, 10, 5, "F(MF,ae)", 0, outfile);
  dpd_buf_close(&F);
  dpd_buf_init(&F, CC_TMP0, 10, 5, 10, 5, 0, "F(MF,ae)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tiaJB", 0, outfile);
  dpd_buf_init(&Z, CC_TMP1, 10, 5, 10, 5, 0, "Z(ib,ae)", 0, outfile);
  dpd_contract222(&T2, &F, &Z, 0, 1, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&F);
  dpd_swap24(&Z, CC_TMP0, 10, 5, "Z(ie,ab)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&Z, CC_TMP0, 10, 5, 10, 5, 0, "Z(ie,ab)", 0, outfile);
  dpd_swap12(&Z, CC_TMP1, 11, 5, "Z(ei,ab)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&Z1, CC_TMP1, 11, 5, 11, 5, 0, "Z(ei,ab)", 0, outfile);
  dpd_swap34(&Z1, CC_TMP0, 11, 5, "Z(ei,ba)", 0, outfile);
  dpd_buf_init(&Z2, CC_TMP0, 11, 5, 11, 5, 0, "Z(ei,ba)", 0, outfile);
  dpd_axpy(&Z2, &Z1, -1.0, 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_init(&W, CC_HBAR, 11, 5, 11, 7, 0, "Weiab", 0, outfile);
  dpd_axpy(&Z1, &W, 1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_close(&W);

  /** -<Mb|Ef> t_Mi^Af - <MA||EF> t_iM^bF + <mA|fE> t_im^bf --> W(Ei,Ab) **/
  dpd_buf_init(&F, CC_FINTS, 10, 5, 10, 5, 0, "F <ia|bc>", 0, outfile);
  dpd_swap24(&F, CC_TMP0, 10, 5, "F(Mf,Eb)", 0, outfile);
  dpd_buf_close(&F);
  dpd_buf_init(&F, CC_TMP0, 10, 5, 10, 5, 0, "F(Mf,Eb)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tIbjA", 0, outfile);
  dpd_buf_init(&Z, CC_TMP1, 10, 5, 10, 5, 0, "Z(iA,Eb)", 0, outfile);
  dpd_contract222(&T2, &F, &Z, 1, 1, -1.0, 0.0, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&F);
  dpd_swap23(&Z, CC_TMP0, 10, 5, "Z(iE,Ab)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&F, CC_FINTS, 10, 5, 10, 5, 0, "F <ia|bc>", 0, outfile);
  dpd_swap23(&F, CC_TMP1, 10, 5, "F(mf,AE)", 0, outfile);
  dpd_buf_close(&F);
  dpd_buf_init(&F, CC_TMP1, 10, 5, 10, 5, 0, "F(mf,AE)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tiajb", 0, outfile);
  dpd_buf_init(&Z, CC_TMP2, 10, 5, 10, 5, 0, "Z(ib,AE)", 0, outfile);
  dpd_contract222(&T2, &F, &Z, 0, 1, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&F);
  dpd_swap24(&Z, CC_TMP1, 10, 5, "Z(iE,Ab)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&F, CC_FINTS, 10, 5, 10, 5, 1, "F <ia|bc>", 0, outfile);
  dpd_swap24(&F, CC_TMP2, 10, 5, "F(MF,EA)", 0, outfile);
  dpd_buf_close(&F);
  dpd_buf_init(&F, CC_TMP2, 10, 5, 10, 5, 0, "F(MF,EA)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tiaJB", 0, outfile); 
  dpd_buf_init(&Z, CC_TMP3, 10, 5, 10, 5, 0, "Z(ib,EA)", 0, outfile);
  dpd_contract222(&T2, &F, &Z, 0, 1, -1.0, 0.0, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&F);
  dpd_swap23(&Z, CC_TMP4, 10, 5, "Z(iE,bA)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&Z, CC_TMP4, 10, 5, 10, 5, 0, "Z(iE,bA)", 0, outfile);
  dpd_swap34(&Z, CC_TMP2, 10, 5, "Z(iE,Ab)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&Z1, CC_TMP0, 10, 5, 10, 5, 0, "Z(iE,Ab)", 0, outfile);
  dpd_buf_init(&Z2, CC_TMP1, 10, 5, 10, 5, 0, "Z(iE,Ab)", 0, outfile);
  dpd_axpy(&Z1, &Z2, 1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_init(&Z1, CC_TMP2, 10, 5, 10, 5, 0, "Z(iE,Ab)", 0, outfile);
  dpd_axpy(&Z1, &Z2, 1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_swap12(&Z2, CC_TMP0, 11, 5, "Z(Ei,Ab)", 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_init(&Z, CC_TMP0, 11, 5, 11, 5, 0, "Z(Ei,Ab)", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 11, 5, 11, 5, 0, "WEiAb", 0, outfile);
  dpd_axpy(&Z, &W, 1.0, 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&W);

  /** -<mB|eF> t_mI^aF - <ma||ef> t_Im^Bf + <Ma|Fe> t_IM^BF --> W(eI,aB) **/
  dpd_buf_init(&F, CC_FINTS, 10, 5, 10, 5, 0, "F <ia|bc>", 0, outfile);
  dpd_swap24(&F, CC_TMP0, 10, 5, "F(mF,eB)", 0, outfile);
  dpd_buf_close(&F);
  dpd_buf_init(&F, CC_TMP0, 10, 5, 10, 5, 0, "F(mF,eB)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tjAIb", 0, outfile);
  dpd_buf_init(&Z, CC_TMP1, 10, 5, 10, 5, 0, "Z(Ia,eB)", 0, outfile);
  dpd_contract222(&T2, &F, &Z, 1, 1, -1.0, 0.0, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&F);
  dpd_swap23(&Z, CC_TMP0, 10, 5, "Z(Ie,aB)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&F, CC_FINTS, 10, 5, 10, 5, 0, "F <ia|bc>", 0, outfile);
  dpd_swap23(&F, CC_TMP1, 10, 5, "F(MF,ae)", 0, outfile);
  dpd_buf_close(&F);
  dpd_buf_init(&F, CC_TMP1, 10, 5, 10, 5, 0, "F(MF,ae)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tIAJB", 0, outfile);
  dpd_buf_init(&Z, CC_TMP2, 10, 5, 10, 5, 0, "Z(IB,ae)", 0, outfile);
  dpd_contract222(&T2, &F, &Z, 0, 1, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&F);
  dpd_swap24(&Z, CC_TMP1, 10, 5, "Z(Ie,aB)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&F, CC_FINTS, 10, 5, 10, 5, 1, "F <ia|bc>", 0, outfile);
  dpd_swap24(&F, CC_TMP2, 10, 5, "F(mf,ea)", 0, outfile);
  dpd_buf_close(&F);
  dpd_buf_init(&F, CC_TMP2, 10, 5, 10, 5, 0, "F(mf,ea)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tIAjb", 0, outfile); 
  dpd_buf_init(&Z, CC_TMP3, 10, 5, 10, 5, 0, "Z(IB,ea)", 0, outfile);
  dpd_contract222(&T2, &F, &Z, 0, 1, -1.0, 0.0, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&F);
  dpd_swap23(&Z, CC_TMP4, 10, 5, "Z(Ie,Ba)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&Z, CC_TMP4, 10, 5, 10, 5, 0, "Z(Ie,Ba)", 0, outfile);
  dpd_swap34(&Z, CC_TMP2, 10, 5, "Z(Ie,aB)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&Z1, CC_TMP0, 10, 5, 10, 5, 0, "Z(Ie,aB)", 0, outfile);
  dpd_buf_init(&Z2, CC_TMP1, 10, 5, 10, 5, 0, "Z(Ie,aB)", 0, outfile);
  dpd_axpy(&Z1, &Z2, 1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_init(&Z1, CC_TMP2, 10, 5, 10, 5, 0, "Z(Ie,aB)", 0, outfile);
  dpd_axpy(&Z1, &Z2, 1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_swap12(&Z2, CC_TMP0, 11, 5, "Z(eI,aB)", 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_init(&Z, CC_TMP0, 11, 5, 11, 5, 0, "Z(eI,aB)", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 11, 5, 11, 5, 0, "WeIaB", 0, outfile);
  dpd_axpy(&Z, &W, 1.0, 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&W);

  /** Final terms of Wabei **/

  /** t_IN^BF <MN||EF> --> Z_IBME **/
  dpd_buf_init(&Z, CC_TMP0, 10, 10, 10, 10, 0, "Z(IB,ME)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tIAJB", 0, outfile);
  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)",
	       0, outfile);
  dpd_contract222(&T2, &D, &Z, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_close(&T2);
  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)",
	       0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tIAjb", 0, outfile);
  dpd_contract222(&T2, &D, &Z, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&D);
  dpd_swap24(&Z, CC_TMP1, 10, 10, "Z(IE,MB)", 0, outfile);
  dpd_buf_close(&Z);

  /** t_M^A ( -<MB||IE> + Z1_IEMB ) --> Z2_EIAB **/
  dpd_buf_init(&Z1, CC_TMP1, 10, 10, 10, 10, 0, "Z(IE,MB)", 0, outfile);
  dpd_buf_init(&C, CC_CINTS, 10, 10, 10, 10, 0, "C <ia||jb>", 0, outfile);
  dpd_axpy(&C, &Z1, -1.0, 0, outfile);
  dpd_buf_close(&C);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_buf_init(&Z2, CC_TMP0, 10, 5, 10, 5, 0, "Z2(IE,AB)", 0, outfile);
  dpd_contract212(&T1, &Z1, &Z2, 0, 2, 1, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&Z1);
  dpd_swap12(&Z2, CC_TMP1, 11, 5, "Z(EI,AB)", 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_init(&Z, CC_TMP1, 11, 5, 11, 5, 0, "Z(EI,AB)", 0, outfile);
  dpd_swap34(&Z, CC_TMP0, 11, 5, "Z(EI,BA)", 0, outfile);
  dpd_buf_close(&Z);

  /** Z1_EIAB - Z2_EIBA --> W_EIAB **/
  dpd_buf_init(&Z1, CC_TMP1, 11, 5, 11, 5, 0, "Z(EI,AB)", 0, outfile);
  dpd_buf_init(&Z2, CC_TMP0, 11, 5, 11, 5, 0, "Z(EI,BA)", 0, outfile);
  dpd_axpy(&Z1, &Z2, -1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_init(&W, CC_HBAR, 11, 5, 11, 7, 0, "WEIAB", 0, outfile);
  dpd_axpy(&Z2, &W, 1.0, 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_close(&W);

  
    /** t_in^bf <mn||ef> --> Z_ibme **/
  dpd_buf_init(&Z, CC_TMP0, 10, 10, 10, 10, 0, "Z(ib,me)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tiajb", 0, outfile);
  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)",
	       0, outfile);
  dpd_contract222(&T2, &D, &Z, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_close(&T2);
  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)",
	       0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tiaJB", 0, outfile);
  dpd_contract222(&T2, &D, &Z, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&D);
  dpd_swap24(&Z, CC_TMP1, 10, 10, "Z(ie,mb)", 0, outfile);
  dpd_buf_close(&Z);

  /** t_m^a ( -<mb||ie> + Z1_iemb ) --> Z2_eiab **/
  dpd_buf_init(&Z1, CC_TMP1, 10, 10, 10, 10, 0, "Z(ie,mb)", 0, outfile);
  dpd_buf_init(&C, CC_CINTS, 10, 10, 10, 10, 0, "C <ia||jb>",
	       0, outfile);
  dpd_axpy(&C, &Z1, -1.0, 0, outfile);
  dpd_buf_close(&C);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_buf_init(&Z2, CC_TMP0, 10, 5, 10, 5, 0, "Z2(ie,ab)", 0, outfile);
  dpd_contract212(&T1, &Z1, &Z2, 0, 2, 1, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&Z1);
  dpd_swap12(&Z2, CC_TMP1, 11, 5, "Z(ei,ab)", 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_init(&Z, CC_TMP1, 11, 5, 11, 5, 0, "Z(ei,ab)", 0, outfile);
  dpd_swap34(&Z, CC_TMP0, 11, 5, "Z(ei,ba)", 0, outfile);
  dpd_buf_close(&Z);

  /** - Z1_eiab + Z2_eiba --> W_eiab **/
  dpd_buf_init(&Z1, CC_TMP1, 11, 5, 11, 5, 0, "Z(ei,ab)", 0, outfile);
  dpd_buf_init(&Z2, CC_TMP0, 11, 5, 11, 5, 0, "Z(ei,ba)", 0, outfile);
  dpd_axpy(&Z1, &Z2, -1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_init(&W, CC_HBAR, 11, 5, 11, 7, 0, "Weiab", 0, outfile);
  dpd_axpy(&Z2, &W, 1.0, 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_close(&W);

  /** t_in^bf  <Mn|Ef> + t_iN^bF <MN||EF> --> Z1_MEib **/
  dpd_buf_init(&Z, CC_TMP0, 10, 10, 10, 10, 0, "Z(ME,ib)", 0, outfile);
  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)",
	       0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tiajb", 0, outfile);
  dpd_contract222(&D, &T2, &Z, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&D);
  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)",
	       0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tiaJB", 0, outfile);
  dpd_contract222(&D, &T2, &Z, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_close(&T2);
  dpd_buf_close(&Z);

  /** - t_Ni^Af <mN|fE> --> Z2_mEiA **/
  dpd_buf_init(&D, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_swap24(&D, CC_TMP1, 10, 11, "D(mE,fN)", 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_init(&D, CC_TMP1, 10, 11, 10, 11, 0, "D(mE,fN)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 0, 5, 0, 5, 0, "tiJaB", 0, outfile);
  dpd_swap24(&T2, CC_TMP2, 10, 11, "T2(iA,fN)", 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_init(&D, CC_TMP1, 10, 11, 10, 11, 0, "D(mE,fN)", 0, outfile);
  dpd_buf_init(&T2, CC_TMP2, 10, 11, 10, 11, 0, "T2(iA,fN)", 0, outfile);
  dpd_buf_init(&Z, CC_TMP3, 10, 10, 10, 10, 0, "Z(mE,iA)", 0, outfile);
  dpd_contract222(&D, &T2, &Z, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&D);
  dpd_buf_close(&Z);

  dpd_buf_init(&Z, CC_TMP0, 10, 10, 10, 10, 0, "Z(ME,ib)", 0, outfile);
  dpd_swap24(&Z, CC_TMP1, 10, 10, "Z(Mb,iE)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&Z, CC_TMP1, 10, 10, 10, 10, 0, "Z(Mb,iE)", 0, outfile);
  dpd_swap34(&Z, CC_TMP0, 10, 11, "Z(Mb,Ei)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&Z, CC_TMP3, 10, 10, 10, 10, 0, "Z(mE,iA)", 0, outfile);
  dpd_swap24(&Z, CC_TMP1, 10, 10, "Z(mA,iE)", 0, outfile);
  dpd_buf_close(&Z);

  /** - t_M^A ( <Mb|Ei> + Z(Mb,Ei) ) --> Z1(Ab,Ei) **/
  dpd_buf_init(&D, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_swap24(&D, CC_TMP2, 10, 11, "D(Mb,Ei)", 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_init(&D, CC_TMP2, 10, 11, 10, 11, 0, "D(Mb,Ei)", 0, outfile);
  dpd_buf_init(&Z, CC_TMP0, 10, 11, 10, 11, 0, "Z(Mb,Ei)", 0, outfile);
  dpd_axpy(&D, &Z, 1.0, 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_init(&Z1, CC_TMP2, 5, 11, 5, 11, 0, "Z1(Ab,Ei)", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_contract212(&T1, &Z, &Z1, 0, 0, 0, -1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&Z);
  dpd_swapbk(&Z1, CC_TMP0, 11, 5, "Z1(Ei,Ab)", 0, outfile);
  dpd_buf_close(&Z1);

  /** t_m^b ( - <mA|iE> + Z(mA,iE) ) --> Z2(Ab,Ei) **/
  dpd_buf_init(&C, CC_CINTS, 10, 10, 10, 10, 0, "C <ia|jb>", 0, outfile);
  dpd_buf_init(&Z, CC_TMP1, 10, 10, 10, 10, 0, "Z(mA,iE)", 0, outfile);
  dpd_axpy(&C, &Z, -1.0, 0, outfile);
  dpd_buf_close(&C);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_buf_init(&Z2, CC_TMP2, 5, 10, 5, 10, 0, "Z2(bA,iE)", 0, outfile);
  dpd_contract212(&T1, &Z, &Z2, 0, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&Z);
  dpd_swap12(&Z2, CC_TMP1, 5, 10, "Z2(Ab,iE)", 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_init(&Z2, CC_TMP1, 5, 10, 5, 10, 0, "Z2(Ab,iE)", 0, outfile);
  dpd_swap34(&Z2, CC_TMP2, 5, 11, "Z2(Ab,Ei)", 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_init(&Z2, CC_TMP2, 5, 11, 5, 11, 0, "Z2(Ab,Ei)", 0, outfile);
  dpd_swapbk(&Z2, CC_TMP1, 11, 5, "Z2(Ei,Ab)", 0, outfile);
  dpd_buf_close(&Z2);

  /** Z1(Ei,Ab) + Z2(Ei,Ab) --> W(Ei,Ab) **/
  dpd_buf_init(&W, CC_HBAR, 11, 5, 11, 5, 0, "WEiAb", 0, outfile);
  dpd_buf_init(&Z1, CC_TMP0, 11, 5, 11, 5, 0, "Z1(Ei,Ab)", 0, outfile);
  dpd_axpy(&Z1, &W, 1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_init(&Z2, CC_TMP1, 11, 5, 11, 5, 0, "Z2(Ei,Ab)",0, outfile);
  dpd_axpy(&Z2, &W, 1.0, 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_close(&W);

  /** t_IN^BF  <mN|eF> + t_In^Bf <mn||ef> --> Z1_meIB **/
  dpd_buf_init(&Z, CC_TMP0, 10, 10, 10, 10, 0, "Z(me,IB)", 0, outfile);
  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)",
	       0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tIAJB", 0, outfile);
  dpd_contract222(&D, &T2, &Z, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&D);
  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)",
	       0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tIAjb", 0, outfile);
  dpd_contract222(&D, &T2, &Z, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_close(&T2);
  dpd_buf_close(&Z);

  /** - t_Ni^Af <mN|fE> --> Z2_mEiA **/
  dpd_buf_init(&D, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_swap24(&D, CC_TMP1, 10, 11, "D(Me,Fn)", 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_init(&D, CC_TMP1, 10, 11, 10, 11, 0, "D(Me,Fn)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 0, 5, 0, 5, 0, "tIjAb", 0, outfile);
  dpd_swap24(&T2, CC_TMP2, 10, 11, "T2(Ia,Fn)", 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_init(&D, CC_TMP1, 10, 11, 10, 11, 0, "D(Me,Fn)", 0, outfile);
  dpd_buf_init(&T2, CC_TMP2, 10, 11, 10, 11, 0, "T2(Ia,Fn)", 0, outfile);
  dpd_buf_init(&Z, CC_TMP3, 10, 10, 10, 10, 0, "Z(Me,Ia)", 0, outfile);
  dpd_contract222(&D, &T2, &Z, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&D);
  dpd_buf_close(&Z);

  dpd_buf_init(&Z, CC_TMP0, 10, 10, 10, 10, 0, "Z(me,IB)", 0, outfile);
  dpd_swap24(&Z, CC_TMP1, 10, 10, "Z(mB,Ie)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&Z, CC_TMP1, 10, 10, 10, 10, 0, "Z(mB,Ie)", 0, outfile);
  dpd_swap34(&Z, CC_TMP0, 10, 11, "Z(mB,eI)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&Z, CC_TMP3, 10, 10, 10, 10, 0, "Z(Me,Ia)", 0, outfile);
  dpd_swap24(&Z, CC_TMP1, 10, 10, "Z(Ma,Ie)", 0, outfile);
  dpd_buf_close(&Z);

  /** - t_m^a ( <mB|eI> + Z(mB,eI) ) --> Z1(aB,eI) **/
  dpd_buf_init(&D, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_swap24(&D, CC_TMP2, 10, 11, "D(mB,eI)", 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_init(&D, CC_TMP2, 10, 11, 10, 11, 0, "D(mB,eI)", 0, outfile);
  dpd_buf_init(&Z, CC_TMP0, 10, 11, 10, 11, 0, "Z(mB,eI)", 0, outfile);
  dpd_axpy(&D, &Z, 1.0, 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_init(&Z1, CC_TMP2, 5, 11, 5, 11, 0, "Z1(aB,eI)", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_contract212(&T1, &Z, &Z1, 0, 0, 0, -1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&Z);
  dpd_swapbk(&Z1, CC_TMP0, 11, 5, "Z1(eI,aB)", 0, outfile);
  dpd_buf_close(&Z1);

  /** t_M^B ( - <Ma|Ie> + Z(Ma,Ie) ) --> Z2(aB,eI) **/
  dpd_buf_init(&C, CC_CINTS, 10, 10, 10, 10, 0, "C <ia|jb>", 0, outfile);
  dpd_buf_init(&Z, CC_TMP1, 10, 10, 10, 10, 0, "Z(Ma,Ie)", 0, outfile);
  dpd_axpy(&C, &Z, -1.0, 0, outfile);
  dpd_buf_close(&C);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_buf_init(&Z2, CC_TMP2, 5, 10, 5, 10, 0, "Z2(Ba,Ie)", 0, outfile);
  dpd_contract212(&T1, &Z, &Z2, 0, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&Z);
  dpd_swap12(&Z2, CC_TMP1, 5, 10, "Z2(aB,Ie)", 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_init(&Z2, CC_TMP1, 5, 10, 5, 10, 0, "Z2(aB,Ie)", 0, outfile);
  dpd_swap34(&Z2, CC_TMP2, 5, 11, "Z2(aB,eI)", 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_init(&Z2, CC_TMP2, 5, 11, 5, 11, 0, "Z2(aB,eI)", 0, outfile);
  dpd_swapbk(&Z2, CC_TMP1, 11, 5, "Z2(eI,aB)", 0, outfile);
  dpd_buf_close(&Z2);

  /** Z1(eI,aB) + Z2(eI,aB) --> W(eI,aB) **/
  dpd_buf_init(&W, CC_HBAR, 11, 5, 11, 5, 0, "WeIaB", 0, outfile);
  dpd_buf_init(&Z1, CC_TMP0, 11, 5, 11, 5, 0, "Z1(eI,aB)", 0, outfile);
  dpd_axpy(&Z1, &W, 1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_init(&Z2, CC_TMP1, 11, 5, 11, 5, 0, "Z2(eI,aB)",0, outfile);
  dpd_axpy(&Z2, &W, 1.0, 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_close(&W);
}
