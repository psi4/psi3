#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void Wmbij_build(void)
{
  struct oe_dpdfile Fme, T1;
  struct dpdbuf W, E, T2, Wmnij, I, Tau, Z, Z1, Z2, C, D;

  dpd_buf_init(&E, CC_EINTS, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)",
	       0, outfile);
  /** <MB||IJ> **/
  dpd_swapbk(&E, CC_HBAR, 10, 2, "WMBIJ", 0, outfile);
  /** <mb||ij> **/
  dpd_swapbk(&E, CC_HBAR, 10, 2, "Wmbij", 0, outfile);
  dpd_buf_close(&E);

  dpd_buf_init(&E, CC_EINTS, 0, 10, 0, 10, 0, "E <ij|ka>", 0, outfile);
  /** <Mb|Ij> **/
  dpd_swapbk(&E, CC_HBAR, 10, 0, "WMbIj", 0, outfile);
  /** <mB|iJ> **/
  dpd_swapbk(&E, CC_HBAR, 10, 0, "WmBiJ", 0, outfile);
  dpd_buf_close(&E);

  /** F_ME t_IJ^EB --> W(MB,IJ) **/
  dpd_oe_file_init(&Fme, CC_OEI, 0, 1, "FME", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 2, 5, 2, 7, 0, "tIJAB", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 2, 10, 2, 0, "WMBIJ", 0, outfile);
  dpd_contract212(&Fme, &T2, &W, 1, 2, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&T2);
  dpd_oe_file_close(&Fme);

  /** F_me t_ij^eb --> W(mb,ij) **/
  dpd_oe_file_init(&Fme, CC_OEI, 0, 1, "Fme", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 2, 5, 2, 7, 0, "tijab", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 2, 10, 2, 0, "Wmbij", 0, outfile);
  dpd_contract212(&Fme, &T2, &W, 1, 2, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&T2);
  dpd_oe_file_close(&Fme);

  /** F_ME t_Ij^Eb --> W(Mb,Ij) **/
  dpd_oe_file_init(&Fme, CC_OEI, 0, 1, "FME", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 0, 5, 0, 5, 0, "tIjAb", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 0, 10, 0, 0, "WMbIj", 0, outfile);
  dpd_contract212(&Fme, &T2, &W, 1, 2, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&T2);
  dpd_oe_file_close(&Fme);

  /** F_me t_iJ^eB --> W(mB,iJ) **/
  dpd_oe_file_init(&Fme, CC_OEI, 0, 1, "Fme", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 0, 5, 0, 5, 0, "tiJaB", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 0, 10, 0, 0, "WmBiJ", 0, outfile);
  dpd_contract212(&Fme, &T2, &W, 1, 2, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&T2);
  dpd_oe_file_close(&Fme);

  /** - t_N^B W_MNIJ --> W(MB,IJ) **/
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_buf_init(&Wmnij, CC_HBAR, 0, 2, 2, 2, 0, "WMNIJ", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 2, 10, 2, 0, "WMBIJ", 0, outfile);
  dpd_contract221(&Wmnij, &T1, &W, 1, 0, 1, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&Wmnij);
  dpd_oe_file_close(&T1);

  /** - t_n^b W_mnij **/
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_buf_init(&Wmnij, CC_HBAR, 0, 2, 2, 2, 0, "Wmnij", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 2, 10, 2, 0, "Wmbij", 0, outfile);
  dpd_contract221(&Wmnij, &T1, &W, 1, 0, 1, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&Wmnij);
  dpd_oe_file_close(&T1);

  /** - t_n^b W_MnIj **/
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_buf_init(&Wmnij, CC_HBAR, 0, 0, 0, 0, 0, "WMnIj", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 0, 10, 0, 0, "WMbIj", 0, outfile);
  dpd_contract221(&Wmnij, &T1, &W, 1, 0, 1, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&Wmnij);
  dpd_oe_file_close(&T1);

  /** - t_N^B W_mNiJ **/
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_buf_init(&Wmnij, CC_HBAR, 0, 0, 0, 0, 0, "WMnIj", 0, outfile);
  dpd_swap12(&Wmnij, CC_TMP0, 0, 0, "WnMIj", 0, outfile);
  dpd_buf_close(&Wmnij);
  dpd_buf_init(&Wmnij, CC_TMP0, 0, 0, 0, 0, 0, "WnMIj", 0, outfile);
  dpd_swap34(&Wmnij, CC_TMP1, 0, 0, "WnMjI", 0, outfile);
  dpd_buf_close(&Wmnij);
  dpd_buf_init(&Wmnij, CC_TMP1, 0, 0, 0, 0, 0, "WnMjI", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 0, 10, 0, 0, "WmBiJ", 0, outfile);
  dpd_contract221(&Wmnij, &T1, &W, 1, 0, 1, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&Wmnij);
  dpd_oe_file_close(&T1);

  /** <MB||EF> tau_IJ^EF **/
  dpd_buf_init(&I, CC_FINTS, 10, 7, 10, 5, 1, "F <ia|bc>", 0, outfile);
  dpd_buf_init(&Tau, CC_TAMPS, 2, 7, 2, 7, 0, "tauIJAB", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 2, 10, 2, 0, "WMBIJ", 0, outfile);
  dpd_contract222(&I, &Tau, &W, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&Tau);
  dpd_buf_close(&I);

  /* <mb||ef> tau_ij^ef **/
  dpd_buf_init(&I, CC_FINTS, 10, 7, 10, 5, 1, "F <ia|bc>", 0, outfile);
  dpd_buf_init(&Tau, CC_TAMPS, 2, 7, 2, 7, 0, "tauijab", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 2, 10, 2, 0, "Wmbij", 0, outfile);
  dpd_contract222(&I, &Tau, &W, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&Tau);
  dpd_buf_close(&I);

  /** <Mb|Ef> tau_Ij^Ef **/
  dpd_buf_init(&I, CC_FINTS, 10, 5, 10, 5, 0, "F <ia|bc>", 0, outfile);
  dpd_buf_init(&Tau, CC_TAMPS, 0, 5, 0, 5, 0, "tauIjAb", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 0, 10, 0, 0, "WMbIj", 0, outfile);
  dpd_contract222(&I, &Tau, &W, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&Tau);
  dpd_buf_close(&I);

  /** <mB|eF> tau_iJ^eF **/
  dpd_buf_init(&I, CC_FINTS, 10, 5, 10, 5, 0, "F <ia|bc>", 0, outfile);
  dpd_buf_init(&Tau, CC_TAMPS, 0, 5, 0, 5, 0, "tauiJaB", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 0, 10, 0, 0, "WmBiJ", 0, outfile);
  dpd_contract222(&I, &Tau, &W, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&Tau);
  dpd_buf_close(&I);

  /* Sort <ij||ka> integrals for the E*T2 contributions */
  dpd_buf_init(&I, CC_EINTS, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)",
	       0, outfile);
  dpd_swap23(&I, CC_TMP0, 0, 10, "I(MI,NE)", 0, outfile);
  dpd_buf_close(&I);

  dpd_buf_init(&I, CC_EINTS, 0, 10, 0, 10, 0, "E <ij|ka>", 0, outfile);
  dpd_swap23(&I, CC_TMP1, 0, 10, "I(MI,NE)", 0, outfile);
  dpd_buf_close(&I);

  /** <MN||IE> t_JN^BE **/
  dpd_buf_init(&I, CC_TMP0, 0, 10, 0, 10, 0, "I(MI,NE)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tIAJB", 0, outfile);
  dpd_buf_init(&Z, CC_TMP2, 0, 10, 0, 10, 0, "Z(MI,JB)", 0, outfile);
  dpd_contract222(&I, &T2, &Z, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&I);

  /** <Mn||Ie> t_Jn^Be **/
  dpd_buf_init(&I, CC_TMP1, 0, 10, 0, 10, 0, "I(MI,NE)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tIAjb", 0, outfile);
  dpd_contract222(&I, &T2, &Z, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&T2);
  dpd_buf_close(&I);

  /** <MN||JE> t_IN^BE **/
  dpd_buf_init(&I, CC_TMP0, 0, 10, 0, 10, 0, "I(MI,NE)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tIAJB", 0, outfile);
  dpd_buf_init(&Z, CC_TMP3, 0, 10, 0, 10, 0, "Z(MJ,IB)", 0, outfile);
  dpd_contract222(&I, &T2, &Z, 0, 0, -1.0, 0.0, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&I);

  /** <Mn||Je> t_In^Be **/
  dpd_buf_init(&I, CC_TMP1, 0, 10, 0, 10, 0, "I(MI,NE)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tIAjb", 0, outfile);
  dpd_contract222(&I, &T2, &Z, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&T2);
  dpd_buf_close(&I);

  dpd_buf_init(&Z1, CC_TMP2, 0, 10, 0, 10, 0, "Z(MI,JB)", 0, outfile);
  dpd_swap23(&Z1, CC_TMP4, 0, 10, "Z(MJ,IB)", 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_init(&Z1, CC_TMP4, 0, 10, 0, 10, 0, "Z(MJ,IB)", 0, outfile);
  dpd_buf_init(&Z2, CC_TMP3, 0, 10, 0, 10, 0, "Z(MJ,IB)", 0, outfile);
  dpd_axpy(&Z1, &Z2, 1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_close(&Z2);
  dpd_buf_init(&Z2, CC_TMP3, 0, 10, 0, 10, 0, "Z(MJ,IB)", 0, outfile);
  dpd_swap24(&Z2, CC_TMP4, 10, 0, "Z(MB,IJ)", 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_init(&Z, CC_TMP4, 10, 0, 10, 0, 0, "Z(MB,IJ)", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 0, 10, 2, 0, "WMBIJ", 0, outfile);
  dpd_axpy(&Z, &W, 1.0, 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&W);

  /** <mn||ie> t_jn^be **/
  dpd_buf_init(&I, CC_TMP0, 0, 10, 0, 10, 0, "I(MI,NE)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tiajb", 0, outfile);
  dpd_buf_init(&Z, CC_TMP2, 0, 10, 0, 10, 0, "Z(mi,jb)", 0, outfile);
  dpd_contract222(&I, &T2, &Z, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&I);

  /** <mN||iE> t_jN^bE **/
  dpd_buf_init(&I, CC_TMP1, 0, 10, 0, 10, 0, "I(MI,NE)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tiaJB", 0, outfile);
  dpd_contract222(&I, &T2, &Z, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&T2);
  dpd_buf_close(&I);

  /** <mn||je> t_in^be **/
  dpd_buf_init(&I, CC_TMP0, 0, 10, 0, 10, 0, "I(MI,NE)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tiajb", 0, outfile);
  dpd_buf_init(&Z, CC_TMP3, 0, 10, 0, 10, 0, "Z(mj,ib)", 0, outfile);
  dpd_contract222(&I, &T2, &Z, 0, 0, -1.0, 0.0, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&I);

  /** <mN||jE> t_iN^bE **/
  dpd_buf_init(&I, CC_TMP1, 0, 10, 0, 10, 0, "I(MI,NE)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tiaJB", 0, outfile);
  dpd_contract222(&I, &T2, &Z, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&T2);
  dpd_buf_close(&I);

  dpd_buf_init(&Z1, CC_TMP2, 0, 10, 0, 10, 0, "Z(mi,jb)", 0, outfile);
  dpd_swap23(&Z1, CC_TMP4, 0, 10, "Z(mj,ib)", 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_init(&Z1, CC_TMP4, 0, 10, 0, 10, 0, "Z(mj,ib)", 0, outfile);
  dpd_buf_init(&Z2, CC_TMP3, 0, 10, 0, 10, 0, "Z(mj,ib)", 0, outfile);
  dpd_axpy(&Z1, &Z2, 1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_close(&Z2);
  dpd_buf_init(&Z2, CC_TMP3, 0, 10, 0, 10, 0, "Z(mj,ib)", 0, outfile);
  dpd_swap24(&Z2, CC_TMP4, 10, 0, "Z(mb,ij)", 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_init(&Z, CC_TMP4, 10, 0, 10, 0, 0, "Z(mb,ij)", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 0, 10, 2, 0, "Wmbij", 0, outfile);
  dpd_axpy(&Z, &W, 1.0, 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&W);

  /** <MN||IE> t_jN^bE **/
  dpd_buf_init(&I, CC_TMP0, 0, 10, 0, 10, 0, "I(MI,NE)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tiaJB", 0, outfile);
  dpd_buf_init(&Z, CC_TMP2, 0, 10, 0, 10, 0, "Z(MI,jb)", 0, outfile);
  dpd_contract222(&I, &T2, &Z, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_swap24(&Z, CC_TMP3, 10, 0, "Z(Mb,jI)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&T2);
  dpd_buf_close(&I);
  dpd_buf_init(&Z, CC_TMP3, 10, 0, 10, 0, 0, "Z(Mb,jI)", 0, outfile);
  dpd_swap34(&Z, CC_TMP4, 10, 0, "Z(Mb,Ij)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&Z, CC_TMP4, 10, 0, 10, 0, 0, "Z(Mb,Ij)", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 0, 10, 0, 0, "WMbIj", 0, outfile);
  dpd_axpy(&Z, &W, 1.0, 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&W);

  /** <Mn|Ie> t_jn^be **/
  dpd_buf_init(&I, CC_TMP1, 0, 10, 0, 10, 0, "I(MI,NE)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tiajb", 0, outfile);
  dpd_buf_init(&Z, CC_TMP2, 0, 10, 0, 10, 0, "Z(MI,jb)", 0, outfile);
  dpd_contract222(&I, &T2, &Z, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_swap24(&Z, CC_TMP3, 10, 0, "Z(Mb,jI)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&T2);
  dpd_buf_close(&I);
  dpd_buf_init(&Z, CC_TMP3, 10, 0, 10, 0, 0, "Z(Mb,jI)", 0, outfile);
  dpd_swap34(&Z, CC_TMP4, 10, 0, "Z(Mb,Ij)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&Z, CC_TMP4, 10, 0, 10, 0, 0, "Z(Mb,Ij)", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 0, 10, 0, 0, "WMbIj", 0, outfile);
  dpd_axpy(&Z, &W, 1.0, 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&W);

  /** <mn||ie> t_Jn^Be **/
  dpd_buf_init(&I, CC_TMP0, 0, 10, 0, 10, 0, "I(MI,NE)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tIAjb", 0, outfile);
  dpd_buf_init(&Z, CC_TMP2, 0, 10, 0, 10, 0, "Z(mi,JB)", 0, outfile);
  dpd_contract222(&I, &T2, &Z, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_swap24(&Z, CC_TMP3, 10, 0, "Z(mB,Ji)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&T2);
  dpd_buf_close(&I);
  dpd_buf_init(&Z, CC_TMP3, 10, 0, 10, 0, 0, "Z(mB,Ji)", 0, outfile);
  dpd_swap34(&Z, CC_TMP4, 10, 0, "Z(mB,iJ)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&Z, CC_TMP4, 10, 0, 10, 0, 0, "Z(mB,iJ)", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 0, 10, 0, 0, "WmBiJ", 0, outfile);
  dpd_axpy(&Z, &W, 1.0, 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&W);

  /** <mN|iE> t_JN^BE **/
  dpd_buf_init(&I, CC_TMP1, 0, 10, 0, 10, 0, "I(MI,NE)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tIAJB", 0, outfile);
  dpd_buf_init(&Z, CC_TMP2, 0, 10, 0, 10, 0, "Z(mi,JB)", 0, outfile);
  dpd_contract222(&I, &T2, &Z, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_swap24(&Z, CC_TMP3, 10, 0, "Z(mB,Ji)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&T2);
  dpd_buf_close(&I);
  dpd_buf_init(&Z, CC_TMP3, 10, 0, 10, 0, 0, "Z(mB,Ji)", 0, outfile);
  dpd_swap34(&Z, CC_TMP4, 10, 0, "Z(mB,iJ)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&Z, CC_TMP4, 10, 0, 10, 0, 0, "Z(mB,iJ)", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 0, 10, 0, 0, "WmBiJ", 0, outfile);
  dpd_axpy(&Z, &W, 1.0, 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&W);

  /* Sort <ai||jk> integrals for remaining E*T2 contributions */
  dpd_buf_init(&I, CC_EINTS, 11, 0, 11, 0, 0, "E <ai|jk>", 0, outfile);
  dpd_swapbk(&I, CC_TMP0, 0, 11, "I(Mn,Ej)", 0, outfile);
  dpd_buf_close(&I);
  dpd_buf_init(&I, CC_TMP0, 0, 11, 0, 11, 0, "I(Mn,Ej)", 0, outfile);
  dpd_swap24(&I, CC_TMP1, 0, 11, "I(Mj,En)", 0, outfile);
  dpd_buf_close(&I);


  /** -<Mn|Ej> t_In^EB **/
  dpd_buf_init(&T2, CC_TAMPS, 0, 5, 0, 5, 0, "tIjAb", 0, outfile);
  dpd_swap24(&T2, CC_TMP0, 10, 11, "T2(Ib,En)", 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_init(&I, CC_TMP1, 0, 11, 0, 11, 0, "I(Mj,En)", 0, outfile);
  dpd_buf_init(&T2, CC_TMP0, 10, 11, 10, 11, 0, "T2(Ib,En)", 0, outfile);
  dpd_buf_init(&Z, CC_TMP2, 0, 10, 0, 10, 0, "Z(Mj,Ib)", 0, outfile);
  dpd_contract222(&I, &T2, &Z, 0, 0, -1.0, 0.0, 0, outfile);
  dpd_swap24(&Z, CC_TMP0, 10, 0, "Z(Mb,Ij)", 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&I);
  dpd_buf_close(&Z);
  dpd_buf_init(&Z, CC_TMP0, 10, 0, 10, 0, 0, "Z(Mb,Ij)", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 0, 10, 0, 0, "WMbIj", 0, outfile);
  dpd_axpy(&Z, &W, 1.0, 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&W);
  

  /** -<mN|eJ> t_iN^eB **/
  dpd_buf_init(&T2, CC_TAMPS, 0, 5, 0, 5, 0, "tiJaB", 0, outfile);
  dpd_swap24(&T2, CC_TMP0, 10, 11, "T2(iB,eN)", 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_init(&I, CC_TMP1, 0, 11, 0, 11, 0, "I(Mj,En)", 0, outfile);
  dpd_buf_init(&T2, CC_TMP0, 10, 11, 10, 11, 0, "T2(iB,eN)", 0, outfile);
  dpd_buf_init(&Z, CC_TMP2, 0, 10, 0, 10, 0, "Z(mJ,iB)", 0, outfile);
  dpd_contract222(&I, &T2, &Z, 0, 0, -1.0, 0.0, 0, outfile);
  dpd_swap24(&Z, CC_TMP0, 10, 0, "Z(mB,iJ)", 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&I);
  dpd_buf_close(&Z);
  dpd_buf_init(&Z, CC_TMP0, 10, 0, 10, 0, 0, "Z(mB,iJ)", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 0, 10, 0, 0, "WmBiJ", 0, outfile);
  dpd_axpy(&Z, &W, 1.0, 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&W);


  /** Prepare intermediates for final term of Wmbij **/

  /** t_JN^BF <MN||EF> --> Z_MBJE **/
  dpd_buf_init(&Z, CC_TMP0, 10, 10, 10, 10, 0, "Z(ME,JB)", 0, outfile);
  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)",
	       0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tIAJB", 0, outfile);
  dpd_contract222(&D, &T2, &Z, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&D);
  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)",
	       0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tIAjb", 0, outfile);
  dpd_contract222(&D, &T2, &Z, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&D);
  dpd_swap24(&Z, CC_TMP1, 10, 10, "Z(MB,JE)", 0, outfile);
  dpd_buf_close(&Z);

  /** t_I^E ( <MB||JE> + Z1_MBJE ) --> Z2_MBIJ **/
  dpd_buf_init(&Z1, CC_TMP1, 10, 10, 10, 10, 0, "Z(MB,JE)", 0, outfile);
  dpd_buf_init(&C, CC_CINTS, 10, 10, 10, 10, 0, "C <ia||jb>", 0, outfile);
  dpd_axpy(&C, &Z1, -1.0, 0, outfile);
  dpd_buf_close(&C);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_buf_init(&Z2, CC_TMP0, 10, 0, 10, 0, 0, "Z1(MB,JI)", 0, outfile);
  dpd_contract221(&Z1, &T1, &Z2, 3, 1, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&Z1);
  dpd_swap34(&Z2, CC_TMP1, 10, 0, "Z2(MB,IJ)", 0, outfile);
  dpd_buf_close(&Z2);

  /** Z1_MBJI(TMP0) - Z2_MBIJ(TMP1) --> W_MBIJ **/
  dpd_buf_init(&Z1, CC_TMP0, 10, 0, 10, 0, 0, "Z1(MB,JI)", 0, outfile);
  dpd_buf_init(&Z2, CC_TMP1, 10, 0, 10, 0, 0, "Z2(MB,IJ)", 0, outfile);
  dpd_axpy(&Z1, &Z2, -1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_init(&W, CC_HBAR, 10, 0, 10, 2, 0, "WMBIJ", 0, outfile);
  dpd_axpy(&Z2, &W, 1.0, 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_close(&W);

  /** t_jn^bf <mn||ef> --> Z_mbje **/
  dpd_buf_init(&Z, CC_TMP0, 10, 10, 10, 10, 0, "Z(me,jb)", 0, outfile);
  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)",
	       0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tiajb", 0, outfile);
  dpd_contract222(&D, &T2, &Z, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&D);
  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)",
	       0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tiaJB", 0, outfile);
  dpd_contract222(&D, &T2, &Z, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&D);
  dpd_swap24(&Z, CC_TMP1, 10, 10, "Z(mb,je)", 0, outfile);
  dpd_buf_close(&Z);

  /** t_i^e ( <mb||je> + Z1_mbje ) --> Z2_mbij **/
  dpd_buf_init(&Z1, CC_TMP1, 10, 10, 10, 10, 0, "Z(mb,je)", 0, outfile);
  dpd_buf_init(&C, CC_CINTS, 10, 10, 10, 10, 0, "C <ia||jb>", 0, outfile);
  dpd_axpy(&C, &Z1, -1.0, 0, outfile);
  dpd_buf_close(&C);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_buf_init(&Z2, CC_TMP0, 10, 0, 10, 0, 0, "Z1(mb,ji)", 0, outfile);
  dpd_contract221(&Z1, &T1, &Z2, 3, 1, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&Z1);
  dpd_swap34(&Z2, CC_TMP1, 10, 0, "Z2(mb,ij)", 0, outfile);
  dpd_buf_close(&Z2);

  /** Z1_mbji(TMP0) - Z2_mbij(TMP1) --> W_mbij **/
  dpd_buf_init(&Z1, CC_TMP0, 10, 0, 10, 0, 0, "Z1(mb,ji)", 0, outfile);
  dpd_buf_init(&Z2, CC_TMP1, 10, 0, 10, 0, 0, "Z2(mb,ij)", 0, outfile);
  dpd_axpy(&Z1, &Z2, -1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_init(&W, CC_HBAR, 10, 0, 10, 2, 0, "Wmbij", 0, outfile);
  dpd_axpy(&Z2, &W, 1.0, 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_close(&W);


  /** <Mn|Ef> t_jn^bf + <MN||EF> t_jN^bF --> Z1_MEjb(TMP0) **/
  dpd_buf_init(&Z1, CC_TMP0, 10, 10, 10, 10, 0, "Z1(ME,jb)", 0, outfile);
  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tiajb", 0, outfile);
  dpd_contract222(&D, &T2, &Z1, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&D);
  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)",
	       0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tiaJB", 0, outfile);
  dpd_contract222(&D, &T2, &Z1, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&D);
  dpd_buf_close(&Z1);

  /** <Mn|Fe> t_In^Fb --> Z2_MeIb (TMP3) **/
  dpd_buf_init(&D, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_swap24(&D, CC_TMP1, 10, 11, "D(Me,Fn)", 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_init(&T2, CC_TAMPS, 0, 5, 0, 5, 0, "tIjAb", 0, outfile);
  dpd_swap24(&T2, CC_TMP2, 10, 11, "T2(Ib,Fn)", 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_init(&Z2, CC_TMP3, 10, 10, 10, 10, 0, "Z(Me,Ib)", 0, outfile);
  dpd_buf_init(&T2, CC_TMP2, 10, 11, 10, 11, 0, "T2(Ib,Fn)", 0, outfile);
  dpd_buf_init(&D, CC_TMP1, 10, 11, 10, 11, 0, "D(Me,Fn)", 0, outfile);
  dpd_contract222(&D, &T2, &Z2, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_close(&T2);
  dpd_buf_close(&Z2);

  dpd_buf_init(&Z1, CC_TMP0, 10, 10, 10, 10, 0, "Z1(ME,jb)", 0, outfile);
  dpd_swap24(&Z1, CC_TMP1, 10, 10, "Z1(Mb,jE)", 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_init(&Z2, CC_TMP3, 10, 10, 10, 10, 0, "Z(Me,Ib)", 0, outfile);
  dpd_swap24(&Z2, CC_TMP0, 10, 10, "Z(Mb,Ie)", 0, outfile);
  dpd_buf_close(&Z2);

  /** t_I^E ( <Mj|Eb> + Z(Mb,jE)(TMP1) ) --> Z(Mb,jI)(TMP1) **/
  dpd_buf_init(&Z1, CC_TMP1, 10, 10, 10, 10, 0, "Z1(Mb,jE)", 0, outfile);
  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)", 0, outfile);
  dpd_axpy(&D, &Z1, 1.0, 0, outfile);
  dpd_buf_close(&D);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_buf_init(&Z, CC_TMP2, 10, 0, 10, 0, 0, "Z(Mb,jI)", 0, outfile);
  dpd_contract221(&Z1, &T1, &Z, 3, 1, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&Z1);
  dpd_swap34(&Z, CC_TMP1, 10, 0, "Z(Mb,Ij)", 0, outfile);
  dpd_buf_close(&Z);

  /** t_j^e ( <Mb|Ie> - Z(Mb,Ie)(TMP0) ) --> Z(Mb,Ij)(TMP2) **/
  dpd_buf_init(&Z1, CC_TMP0, 10, 10, 10, 10, 0, "Z(Mb,Ie)", 0, outfile);
  dpd_buf_init(&C, CC_CINTS, 10, 10, 10, 10, 0, "C <ia|jb>", 0, outfile);
  dpd_axpy(&C, &Z1, -1.0, 0, outfile);
  dpd_buf_close(&C);
  dpd_scm(&Z1, -1.0, 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_buf_init(&Z, CC_TMP2, 10, 0, 10, 0, 0, "Z(Mb,Ij)", 0, outfile);
  dpd_contract221(&Z1, &T1, &Z, 3, 1, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&Z1);
  dpd_buf_close(&Z);

  /** Z(Mb,Ij) (TMP1) + Z(Mb,Ij) (TMP2) --> W_MbIj **/
  dpd_buf_init(&Z1, CC_TMP1, 10, 0, 10, 0, 0, "Z(Mb,Ij)", 0, outfile);
  dpd_buf_init(&Z2, CC_TMP2, 10, 0, 10, 0, 0, "Z(Mb,Ij)", 0, outfile);
  dpd_axpy(&Z2, &Z1, 1.0, 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_init(&W, CC_HBAR, 10, 0, 10, 0, 0, "WMbIj", 0, outfile);
  dpd_axpy(&Z1, &W, 1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_close(&W);

  /** t_JN^BF <mN|eF> + t_Jn^Bf <mn||ef> --> Z(mB,Je) (TMP1) **/
  dpd_buf_init(&Z, CC_TMP0, 10, 10, 10, 10, 0, "Z(me,JB)", 0, outfile);
  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)",
	       0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tIAJB", 0, outfile);
  dpd_contract222(&D, &T2, &Z, 0 , 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&D);
  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)",
	       0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tIAjb", 0, outfile);
  dpd_contract222(&D, &T2, &Z, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&D);
  dpd_swap24(&Z, CC_TMP1, 10, 10, "Z(mB,Je)", 0, outfile);
  dpd_buf_close(&Z);

  /** -t_Ni^Bf <mN|fE> --> Z(mB,iE) (TMP0) **/
  dpd_buf_init(&Z, CC_TMP2, 10, 10, 10, 10, 0, "Z(mE,iB)", 0, outfile);
  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)",
	       0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 10, 10, 10, 10, 0, "tjAIb", 0, outfile);
  dpd_contract222(&D, &T2, &Z, 0, 0, -1.0, 0.0, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&D);
  dpd_swap24(&Z, CC_TMP0, 10, 10, "Z(mB,iE)", 0, outfile);
  dpd_buf_close(&Z);

  /** t_i^e ( <mJ|eB> + Z(mB,Je) ) --> Z1(mB,iJ) (TMP1) **/
  dpd_buf_init(&Z, CC_TMP1, 10, 10, 10, 10, 0, "Z(mB,Je)", 0, outfile);
  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)", 0, outfile);
  dpd_axpy(&D, &Z, 1.0, 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_init(&Z1, CC_TMP2, 10, 0, 10, 0, 0, "Z(mB,Ji)", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_contract221(&Z, &T1, &Z1, 3, 1, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&Z);
  dpd_swap34(&Z1, CC_TMP1, 10, 0, "Z1(mB,iJ)", 0, outfile);
  dpd_buf_close(&Z1);

  /** t_J^E ( <mB|iE> + Z(mB,iE) ) + Z1(mB,Ij) (TMP1) --> Z2(mB,iJ) (TMP2) **/
  dpd_buf_init(&Z, CC_TMP0, 10, 10, 10, 10, 0, "Z(mB,iE)", 0, outfile);
  dpd_buf_init(&C, CC_CINTS, 10, 10, 10, 10, 0, "C <ia|jb>", 0, outfile);
  dpd_axpy(&C, &Z, 1.0, 0, outfile);
  dpd_buf_close(&C);
  dpd_buf_init(&Z2, CC_TMP2, 10, 0, 10, 0, 0, "Z2(mB,iJ)", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_contract221(&Z, &T1, &Z2, 3, 1, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&Z);
  dpd_buf_init(&Z1, CC_TMP1, 10, 0, 10, 0, 0, "Z1(mB,iJ)", 0, outfile);
  dpd_axpy(&Z1, &Z2, 1.0, 0, outfile);
  dpd_buf_close(&Z1);

  /** Z2(mB,iJ) --> W_mBiJ **/
  dpd_buf_init(&W, CC_HBAR, 10, 0, 10, 0, 0, "WmBiJ", 0, outfile);
  dpd_axpy(&Z2, &W, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&Z2);
}
