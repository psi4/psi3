#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void Wmbij_build(void)
{
  dpdfile2 Fme, T1;
  dpdbuf4 W, E, T2, Wmnij, I, Tau, Z, Z1, Z2, C, D;

  dpd_buf4_init(&E, CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
  /** <MB||IJ> **/
  dpd_buf4_sort(&E, CC_HBAR, rspq, 10, 2, "WMBIJ");
  /** <mb||ij> **/
  dpd_buf4_sort(&E, CC_HBAR, rspq, 10, 2, "Wmbij");
  dpd_buf4_close(&E);

  dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  /** <Mb|Ij> **/
  dpd_buf4_sort(&E, CC_HBAR, rspq, 10, 0, "WMbIj");
  /** <mB|iJ> **/
  dpd_buf4_sort(&E, CC_HBAR, rspq, 10, 0, "WmBiJ");
  dpd_buf4_close(&E);

  /** F_ME t_IJ^EB --> W(MB,IJ) **/
  dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "FME");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 2, 10, 2, 0, "WMBIJ");
  dpd_contract244(&Fme, &T2, &W, 1, 2, 0, 1.0, 1.0);
  dpd_buf4_close(&W);
  dpd_buf4_close(&T2);
  dpd_file2_close(&Fme);

  /** F_me t_ij^eb --> W(mb,ij) **/
  dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "Fme");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 2, 10, 2, 0, "Wmbij");
  dpd_contract244(&Fme, &T2, &W, 1, 2, 0, 1.0, 1.0);
  dpd_buf4_close(&W);
  dpd_buf4_close(&T2);
  dpd_file2_close(&Fme);

  /** F_ME t_Ij^Eb --> W(Mb,Ij) **/
  dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "FME");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
  dpd_contract244(&Fme, &T2, &W, 1, 2, 0, 1.0, 1.0);
  dpd_buf4_close(&W);
  dpd_buf4_close(&T2);
  dpd_file2_close(&Fme);

  /** F_me t_iJ^eB --> W(mB,iJ) **/
  dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "Fme");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 0, 0, "WmBiJ");
  dpd_contract244(&Fme, &T2, &W, 1, 2, 0, 1.0, 1.0);
  dpd_buf4_close(&W);
  dpd_buf4_close(&T2);
  dpd_file2_close(&Fme);

  /** - t_N^B W_MNIJ --> W(MB,IJ) **/
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_buf4_init(&Wmnij, CC_HBAR, 0, 0, 2, 2, 2, 0, "WMNIJ");
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 2, 10, 2, 0, "WMBIJ");
  dpd_contract424(&Wmnij, &T1, &W, 1, 0, 1, -1.0, 1.0);
  dpd_buf4_close(&W);
  dpd_buf4_close(&Wmnij);
  dpd_file2_close(&T1);

  /** - t_n^b W_mnij **/
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_buf4_init(&Wmnij, CC_HBAR, 0, 0, 2, 2, 2, 0, "Wmnij");
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 2, 10, 2, 0, "Wmbij");
  dpd_contract424(&Wmnij, &T1, &W, 1, 0, 1, -1.0, 1.0);
  dpd_buf4_close(&W);
  dpd_buf4_close(&Wmnij);
  dpd_file2_close(&T1);

  /** - t_n^b W_MnIj **/
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_buf4_init(&Wmnij, CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
  dpd_contract424(&Wmnij, &T1, &W, 1, 0, 1, -1.0, 1.0);
  dpd_buf4_close(&W);
  dpd_buf4_close(&Wmnij);
  dpd_file2_close(&T1);

  /** - t_N^B W_mNiJ **/
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_buf4_init(&Wmnij, CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
  dpd_buf4_sort(&Wmnij, CC_TMP0, qprs, 0, 0, "WnMIj");
  dpd_buf4_close(&Wmnij);
  dpd_buf4_init(&Wmnij, CC_TMP0, 0, 0, 0, 0, 0, 0, "WnMIj");
  dpd_buf4_sort(&Wmnij, CC_TMP1, pqsr, 0, 0, "WnMjI");
  dpd_buf4_close(&Wmnij);
  dpd_buf4_init(&Wmnij, CC_TMP1, 0, 0, 0, 0, 0, 0, "WnMjI");
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 0, 0, "WmBiJ");
  dpd_contract424(&Wmnij, &T1, &W, 1, 0, 1, -1.0, 1.0);
  dpd_buf4_close(&W);
  dpd_buf4_close(&Wmnij);
  dpd_file2_close(&T1);

  /** <MB||EF> tau_IJ^EF **/
  dpd_buf4_init(&I, CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
  dpd_buf4_init(&Tau, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 2, 10, 2, 0, "WMBIJ");
  dpd_contract444(&I, &Tau, &W, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&W);
  dpd_buf4_close(&Tau);
  dpd_buf4_close(&I);

  /* <mb||ef> tau_ij^ef **/
  dpd_buf4_init(&I, CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
  dpd_buf4_init(&Tau, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 2, 10, 2, 0, "Wmbij");
  dpd_contract444(&I, &Tau, &W, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&W);
  dpd_buf4_close(&Tau);
  dpd_buf4_close(&I);

  /** <Mb|Ef> tau_Ij^Ef **/
  dpd_buf4_init(&I, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_buf4_init(&Tau, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
  dpd_contract444(&I, &Tau, &W, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&W);
  dpd_buf4_close(&Tau);
  dpd_buf4_close(&I);

  /** <mB|eF> tau_iJ^eF **/
  dpd_buf4_init(&I, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_buf4_init(&Tau, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauiJaB");
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 0, 0, "WmBiJ");
  dpd_contract444(&I, &Tau, &W, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&W);
  dpd_buf4_close(&Tau);
  dpd_buf4_close(&I);

  /* Sort <ij||ka> integrals for the E*T2 contributions */
  dpd_buf4_init(&I, CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
  dpd_buf4_sort(&I, CC_TMP0, prqs, 0, 10, "I(MI,NE)");
  dpd_buf4_close(&I);

  dpd_buf4_init(&I, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  dpd_buf4_sort(&I, CC_TMP1, prqs, 0, 10, "I(MI,NE)");
  dpd_buf4_close(&I);

  /** <MN||IE> t_JN^BE **/
  dpd_buf4_init(&I, CC_TMP0, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
  dpd_buf4_init(&Z, CC_TMP2, 0, 0, 10, 0, 10, 0, "Z(MI,JB)");
  dpd_contract444(&I, &T2, &Z, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&I);

  /** <Mn||Ie> t_Jn^Be **/
  dpd_buf4_init(&I, CC_TMP1, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
  dpd_contract444(&I, &T2, &Z, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&I);

  /** <MN||JE> t_IN^BE **/
  dpd_buf4_init(&I, CC_TMP0, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
  dpd_buf4_init(&Z, CC_TMP3, 0, 0, 10, 0, 10, 0, "Z(MJ,IB)");
  dpd_contract444(&I, &T2, &Z, 0, 0, -1.0, 0.0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&I);

  /** <Mn||Je> t_In^Be **/
  dpd_buf4_init(&I, CC_TMP1, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
  dpd_contract444(&I, &T2, &Z, 0, 0, -1.0, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&I);

  dpd_buf4_init(&Z1, CC_TMP2, 0, 0, 10, 0, 10, 0, "Z(MI,JB)");
  dpd_buf4_sort(&Z1, CC_TMP4, prqs, 0, 10, "Z(MJ,IB)");
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z1, CC_TMP4, 0, 0, 10, 0, 10, 0, "Z(MJ,IB)");
  dpd_buf4_init(&Z2, CC_TMP3, 0, 0, 10, 0, 10, 0, "Z(MJ,IB)");
  dpd_buf4_axpy(&Z1, &Z2, 1.0);
  dpd_buf4_close(&Z1);
  dpd_buf4_close(&Z2);
  dpd_buf4_init(&Z2, CC_TMP3, 0, 0, 10, 0, 10, 0, "Z(MJ,IB)");
  dpd_buf4_sort(&Z2, CC_TMP4, psrq, 10, 0, "Z(MB,IJ)");
  dpd_buf4_close(&Z2);
  dpd_buf4_init(&Z, CC_TMP4, 0, 10, 0, 10, 0, 0, "Z(MB,IJ)");
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 2, 0, "WMBIJ");
  dpd_buf4_axpy(&Z, &W, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&W);

  /** <mn||ie> t_jn^be **/
  dpd_buf4_init(&I, CC_TMP0, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
  dpd_buf4_init(&Z, CC_TMP2, 0, 0, 10, 0, 10, 0, "Z(mi,jb)");
  dpd_contract444(&I, &T2, &Z, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&I);

  /** <mN||iE> t_jN^bE **/
  dpd_buf4_init(&I, CC_TMP1, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
  dpd_contract444(&I, &T2, &Z, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&I);

  /** <mn||je> t_in^be **/
  dpd_buf4_init(&I, CC_TMP0, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
  dpd_buf4_init(&Z, CC_TMP3, 0, 0, 10, 0, 10, 0, "Z(mj,ib)");
  dpd_contract444(&I, &T2, &Z, 0, 0, -1.0, 0.0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&I);

  /** <mN||jE> t_iN^bE **/
  dpd_buf4_init(&I, CC_TMP1, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
  dpd_contract444(&I, &T2, &Z, 0, 0, -1.0, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&I);

  dpd_buf4_init(&Z1, CC_TMP2, 0, 0, 10, 0, 10, 0, "Z(mi,jb)");
  dpd_buf4_sort(&Z1, CC_TMP4, prqs, 0, 10, "Z(mj,ib)");
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z1, CC_TMP4, 0, 0, 10, 0, 10, 0, "Z(mj,ib)");
  dpd_buf4_init(&Z2, CC_TMP3, 0, 0, 10, 0, 10, 0, "Z(mj,ib)");
  dpd_buf4_axpy(&Z1, &Z2, 1.0);
  dpd_buf4_close(&Z1);
  dpd_buf4_close(&Z2);
  dpd_buf4_init(&Z2, CC_TMP3, 0, 0, 10, 0, 10, 0, "Z(mj,ib)");
  dpd_buf4_sort(&Z2, CC_TMP4, psrq, 10, 0, "Z(mb,ij)");
  dpd_buf4_close(&Z2);
  dpd_buf4_init(&Z, CC_TMP4, 0, 10, 0, 10, 0, 0, "Z(mb,ij)");
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 2, 0, "Wmbij");
  dpd_buf4_axpy(&Z, &W, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&W);

  /** <MN||IE> t_jN^bE **/
  dpd_buf4_init(&I, CC_TMP0, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
  dpd_buf4_init(&Z, CC_TMP2, 0, 0, 10, 0, 10, 0, "Z(MI,jb)");
  dpd_contract444(&I, &T2, &Z, 0, 0, 1.0, 0.0);
  dpd_buf4_sort(&Z, CC_TMP3, psrq, 10, 0, "Z(Mb,jI)");
  dpd_buf4_close(&Z);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&I);
  dpd_buf4_init(&Z, CC_TMP3, 0, 10, 0, 10, 0, 0, "Z(Mb,jI)");
  dpd_buf4_sort(&Z, CC_TMP4, pqsr, 10, 0, "Z(Mb,Ij)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, CC_TMP4, 0, 10, 0, 10, 0, 0, "Z(Mb,Ij)");
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
  dpd_buf4_axpy(&Z, &W, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&W);

  /** <Mn|Ie> t_jn^be **/
  dpd_buf4_init(&I, CC_TMP1, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
  dpd_buf4_init(&Z, CC_TMP2, 0, 0, 10, 0, 10, 0, "Z(MI,jb)");
  dpd_contract444(&I, &T2, &Z, 0, 0, 1.0, 0.0);
  dpd_buf4_sort(&Z, CC_TMP3, psrq, 10, 0, "Z(Mb,jI)");
  dpd_buf4_close(&Z);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&I);
  dpd_buf4_init(&Z, CC_TMP3, 0, 10, 0, 10, 0, 0, "Z(Mb,jI)");
  dpd_buf4_sort(&Z, CC_TMP4, pqsr, 10, 0, "Z(Mb,Ij)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, CC_TMP4, 0, 10, 0, 10, 0, 0, "Z(Mb,Ij)");
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
  dpd_buf4_axpy(&Z, &W, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&W);

  /** <mn||ie> t_Jn^Be **/
  dpd_buf4_init(&I, CC_TMP0, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
  dpd_buf4_init(&Z, CC_TMP2, 0, 0, 10, 0, 10, 0, "Z(mi,JB)");
  dpd_contract444(&I, &T2, &Z, 0, 0, 1.0, 0.0);
  dpd_buf4_sort(&Z, CC_TMP3, psrq, 10, 0, "Z(mB,Ji)");
  dpd_buf4_close(&Z);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&I);
  dpd_buf4_init(&Z, CC_TMP3, 0, 10, 0, 10, 0, 0, "Z(mB,Ji)");
  dpd_buf4_sort(&Z, CC_TMP4, pqsr, 10, 0, "Z(mB,iJ)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, CC_TMP4, 0, 10, 0, 10, 0, 0, "Z(mB,iJ)");
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 0, 0, "WmBiJ");
  dpd_buf4_axpy(&Z, &W, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&W);

  /** <mN|iE> t_JN^BE **/
  dpd_buf4_init(&I, CC_TMP1, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
  dpd_buf4_init(&Z, CC_TMP2, 0, 0, 10, 0, 10, 0, "Z(mi,JB)");
  dpd_contract444(&I, &T2, &Z, 0, 0, 1.0, 0.0);
  dpd_buf4_sort(&Z, CC_TMP3, psrq, 10, 0, "Z(mB,Ji)");
  dpd_buf4_close(&Z);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&I);
  dpd_buf4_init(&Z, CC_TMP3, 0, 10, 0, 10, 0, 0, "Z(mB,Ji)");
  dpd_buf4_sort(&Z, CC_TMP4, pqsr, 10, 0, "Z(mB,iJ)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, CC_TMP4, 0, 10, 0, 10, 0, 0, "Z(mB,iJ)");
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 0, 0, "WmBiJ");
  dpd_buf4_axpy(&Z, &W, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&W);

  /* Sort <ai||jk> integrals for remaining E*T2 contributions */
  dpd_buf4_init(&I, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
  dpd_buf4_sort(&I, CC_TMP0, rspq, 0, 11, "I(Mn,Ej)");
  dpd_buf4_close(&I);
  dpd_buf4_init(&I, CC_TMP0, 0, 0, 11, 0, 11, 0, "I(Mn,Ej)");
  dpd_buf4_sort(&I, CC_TMP1, psrq, 0, 11, "I(Mj,En)");
  dpd_buf4_close(&I);


  /** -<Mn|Ej> t_In^EB **/
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_sort(&T2, CC_TMP0, psrq, 10, 11, "T2(Ib,En)");
  dpd_buf4_close(&T2);
  dpd_buf4_init(&I, CC_TMP1, 0, 0, 11, 0, 11, 0, "I(Mj,En)");
  dpd_buf4_init(&T2, CC_TMP0, 0, 10, 11, 10, 11, 0, "T2(Ib,En)");
  dpd_buf4_init(&Z, CC_TMP2, 0, 0, 10, 0, 10, 0, "Z(Mj,Ib)");
  dpd_contract444(&I, &T2, &Z, 0, 0, -1.0, 0.0);
  dpd_buf4_sort(&Z, CC_TMP0, psrq, 10, 0, "Z(Mb,Ij)");
  dpd_buf4_close(&T2);
  dpd_buf4_close(&I);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, CC_TMP0, 0, 10, 0, 10, 0, 0, "Z(Mb,Ij)");
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
  dpd_buf4_axpy(&Z, &W, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&W);
  

  /** -<mN|eJ> t_iN^eB **/
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
  dpd_buf4_sort(&T2, CC_TMP0, psrq, 10, 11, "T2(iB,eN)");
  dpd_buf4_close(&T2);
  dpd_buf4_init(&I, CC_TMP1, 0, 0, 11, 0, 11, 0, "I(Mj,En)");
  dpd_buf4_init(&T2, CC_TMP0, 0, 10, 11, 10, 11, 0, "T2(iB,eN)");
  dpd_buf4_init(&Z, CC_TMP2, 0, 0, 10, 0, 10, 0, "Z(mJ,iB)");
  dpd_contract444(&I, &T2, &Z, 0, 0, -1.0, 0.0);
  dpd_buf4_sort(&Z, CC_TMP0, psrq, 10, 0, "Z(mB,iJ)");
  dpd_buf4_close(&T2);
  dpd_buf4_close(&I);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, CC_TMP0, 0, 10, 0, 10, 0, 0, "Z(mB,iJ)");
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 0, 0, "WmBiJ");
  dpd_buf4_axpy(&Z, &W, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&W);


  /** Prepare intermediates for final term of Wmbij **/

  /** t_JN^BF <MN||EF> --> Z_MBJE **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(ME,JB)");
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
  dpd_contract444(&D, &T2, &Z, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&D);
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
  dpd_contract444(&D, &T2, &Z, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&D);
  dpd_buf4_sort(&Z, CC_TMP1, psrq, 10, 10, "Z(MB,JE)");
  dpd_buf4_close(&Z);

  /** t_I^E ( <MB||JE> + Z1_MBJE ) --> Z2_MBIJ **/
  dpd_buf4_init(&Z1, CC_TMP1, 0, 10, 10, 10, 10, 0, "Z(MB,JE)");
  dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
  dpd_buf4_axpy(&C, &Z1, -1.0);
  dpd_buf4_close(&C);
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_buf4_init(&Z2, CC_TMP0, 0, 10, 0, 10, 0, 0, "Z1(MB,JI)");
  dpd_contract424(&Z1, &T1, &Z2, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z1);
  dpd_buf4_sort(&Z2, CC_TMP1, pqsr, 10, 0, "Z2(MB,IJ)");
  dpd_buf4_close(&Z2);

  /** Z1_MBJI(TMP0) - Z2_MBIJ(TMP1) --> W_MBIJ **/
  dpd_buf4_init(&Z1, CC_TMP0, 0, 10, 0, 10, 0, 0, "Z1(MB,JI)");
  dpd_buf4_init(&Z2, CC_TMP1, 0, 10, 0, 10, 0, 0, "Z2(MB,IJ)");
  dpd_buf4_axpy(&Z1, &Z2, -1.0);
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 2, 0, "WMBIJ");
  dpd_buf4_axpy(&Z2, &W, 1.0);
  dpd_buf4_close(&Z2);
  dpd_buf4_close(&W);

  /** t_jn^bf <mn||ef> --> Z_mbje **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(me,jb)");
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
  dpd_contract444(&D, &T2, &Z, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&D);
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
  dpd_contract444(&D, &T2, &Z, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&D);
  dpd_buf4_sort(&Z, CC_TMP1, psrq, 10, 10, "Z(mb,je)");
  dpd_buf4_close(&Z);

  /** t_i^e ( <mb||je> + Z1_mbje ) --> Z2_mbij **/
  dpd_buf4_init(&Z1, CC_TMP1, 0, 10, 10, 10, 10, 0, "Z(mb,je)");
  dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
  dpd_buf4_axpy(&C, &Z1, -1.0);
  dpd_buf4_close(&C);
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_buf4_init(&Z2, CC_TMP0, 0, 10, 0, 10, 0, 0, "Z1(mb,ji)");
  dpd_contract424(&Z1, &T1, &Z2, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z1);
  dpd_buf4_sort(&Z2, CC_TMP1, pqsr, 10, 0, "Z2(mb,ij)");
  dpd_buf4_close(&Z2);

  /** Z1_mbji(TMP0) - Z2_mbij(TMP1) --> W_mbij **/
  dpd_buf4_init(&Z1, CC_TMP0, 0, 10, 0, 10, 0, 0, "Z1(mb,ji)");
  dpd_buf4_init(&Z2, CC_TMP1, 0, 10, 0, 10, 0, 0, "Z2(mb,ij)");
  dpd_buf4_axpy(&Z1, &Z2, -1.0);
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 2, 0, "Wmbij");
  dpd_buf4_axpy(&Z2, &W, 1.0);
  dpd_buf4_close(&Z2);
  dpd_buf4_close(&W);


  /** <Mn|Ef> t_jn^bf + <MN||EF> t_jN^bF --> Z1_MEjb(TMP0) **/
  dpd_buf4_init(&Z1, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z1(ME,jb)");
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
  dpd_contract444(&D, &T2, &Z1, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&D);
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
  dpd_contract444(&D, &T2, &Z1, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&D);
  dpd_buf4_close(&Z1);

  /** <Mn|Fe> t_In^Fb --> Z2_MeIb (TMP3) **/
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_sort(&D, CC_TMP1, psrq, 10, 11, "D(Me,Fn)");
  dpd_buf4_close(&D);
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_sort(&T2, CC_TMP2, psrq, 10, 11, "T2(Ib,Fn)");
  dpd_buf4_close(&T2);
  dpd_buf4_init(&Z2, CC_TMP3, 0, 10, 10, 10, 10, 0, "Z(Me,Ib)");
  dpd_buf4_init(&T2, CC_TMP2, 0, 10, 11, 10, 11, 0, "T2(Ib,Fn)");
  dpd_buf4_init(&D, CC_TMP1, 0, 10, 11, 10, 11, 0, "D(Me,Fn)");
  dpd_contract444(&D, &T2, &Z2, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&D);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&Z2);

  dpd_buf4_init(&Z1, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z1(ME,jb)");
  dpd_buf4_sort(&Z1, CC_TMP1, psrq, 10, 10, "Z1(Mb,jE)");
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z2, CC_TMP3, 0, 10, 10, 10, 10, 0, "Z(Me,Ib)");
  dpd_buf4_sort(&Z2, CC_TMP0, psrq, 10, 10, "Z(Mb,Ie)");
  dpd_buf4_close(&Z2);

  /** t_I^E ( <Mj|Eb> + Z(Mb,jE)(TMP1) ) --> Z(Mb,jI)(TMP1) **/
  dpd_buf4_init(&Z1, CC_TMP1, 0, 10, 10, 10, 10, 0, "Z1(Mb,jE)");
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
  dpd_buf4_axpy(&D, &Z1, 1.0);
  dpd_buf4_close(&D);
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_buf4_init(&Z, CC_TMP2, 0, 10, 0, 10, 0, 0, "Z(Mb,jI)");
  dpd_contract424(&Z1, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z1);
  dpd_buf4_sort(&Z, CC_TMP1, pqsr, 10, 0, "Z(Mb,Ij)");
  dpd_buf4_close(&Z);

  /** t_j^e ( <Mb|Ie> - Z(Mb,Ie)(TMP0) ) --> Z(Mb,Ij)(TMP2) **/
  dpd_buf4_init(&Z1, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(Mb,Ie)");
  dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
  dpd_buf4_axpy(&C, &Z1, -1.0);
  dpd_buf4_close(&C);
  dpd_buf4_scm(&Z1, -1.0);
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_buf4_init(&Z, CC_TMP2, 0, 10, 0, 10, 0, 0, "Z(Mb,Ij)");
  dpd_contract424(&Z1, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z1);
  dpd_buf4_close(&Z);

  /** Z(Mb,Ij) (TMP1) + Z(Mb,Ij) (TMP2) --> W_MbIj **/
  dpd_buf4_init(&Z1, CC_TMP1, 0, 10, 0, 10, 0, 0, "Z(Mb,Ij)");
  dpd_buf4_init(&Z2, CC_TMP2, 0, 10, 0, 10, 0, 0, "Z(Mb,Ij)");
  dpd_buf4_axpy(&Z2, &Z1, 1.0);
  dpd_buf4_close(&Z2);
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
  dpd_buf4_axpy(&Z1, &W, 1.0);
  dpd_buf4_close(&Z1);
  dpd_buf4_close(&W);

  /** t_JN^BF <mN|eF> + t_Jn^Bf <mn||ef> --> Z(mB,Je) (TMP1) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(me,JB)");
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
  dpd_contract444(&D, &T2, &Z, 0 , 0, 1.0, 0.0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&D);
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
  dpd_contract444(&D, &T2, &Z, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&D);
  dpd_buf4_sort(&Z, CC_TMP1, psrq, 10, 10, "Z(mB,Je)");
  dpd_buf4_close(&Z);

  /** -t_Ni^Bf <mN|fE> --> Z(mB,iE) (TMP0) **/
  dpd_buf4_init(&Z, CC_TMP2, 0, 10, 10, 10, 10, 0, "Z(mE,iB)");
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tjAIb");
  dpd_contract444(&D, &T2, &Z, 0, 0, -1.0, 0.0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&D);
  dpd_buf4_sort(&Z, CC_TMP0, psrq, 10, 10, "Z(mB,iE)");
  dpd_buf4_close(&Z);

  /** t_i^e ( <mJ|eB> + Z(mB,Je) ) --> Z1(mB,iJ) (TMP1) **/
  dpd_buf4_init(&Z, CC_TMP1, 0, 10, 10, 10, 10, 0, "Z(mB,Je)");
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
  dpd_buf4_axpy(&D, &Z, 1.0);
  dpd_buf4_close(&D);
  dpd_buf4_init(&Z1, CC_TMP2, 0, 10, 0, 10, 0, 0, "Z(mB,Ji)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&Z, &T1, &Z1, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_sort(&Z1, CC_TMP1, pqsr, 10, 0, "Z1(mB,iJ)");
  dpd_buf4_close(&Z1);

  /** t_J^E ( <mB|iE> + Z(mB,iE) ) + Z1(mB,Ij) (TMP1) --> Z2(mB,iJ) (TMP2) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(mB,iE)");
  dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
  dpd_buf4_axpy(&C, &Z, 1.0);
  dpd_buf4_close(&C);
  dpd_buf4_init(&Z2, CC_TMP2, 0, 10, 0, 10, 0, 0, "Z2(mB,iJ)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&Z, &T1, &Z2, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z1, CC_TMP1, 0, 10, 0, 10, 0, 0, "Z1(mB,iJ)");
  dpd_buf4_axpy(&Z1, &Z2, 1.0);
  dpd_buf4_close(&Z1);

  /** Z2(mB,iJ) --> W_mBiJ **/
  dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 0, 0, "WmBiJ");
  dpd_buf4_axpy(&Z2, &W, 1.0);
  dpd_buf4_close(&W);
  dpd_buf4_close(&Z2);
}
