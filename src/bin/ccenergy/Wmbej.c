#include <stdio.h>
#include <dpd.h>
#include <qt.h>
#define EXTERN
#include "globals.h"

void Wmbej_build(void)
{
  struct dpdbuf WMBEJ, Wmbej, WMbEj, WmBeJ, WmBEj, WMbeJ; 
  struct dpdbuf C, D, E, F, X, tIAjb, tiaJB;
  struct oe_dpdfile tIA, tia;

/*  timer_on("Wmbej_build", outfile); */

  dpd_buf_init(&C, CC_CINTS, 10, 10, 10, 10, 0, "C <ia||jb>", 0, outfile);
  dpd_copy(&C, CC_HBAR, "WMBEJ", 0, outfile);
  dpd_buf_close(&C);
  dpd_buf_init(&WMBEJ,CC_HBAR, 10, 10, 10, 10, 0, "WMBEJ", 0, outfile);
  dpd_scm(&WMBEJ, -1, 0, outfile);
  dpd_buf_close(&WMBEJ);

  dpd_buf_init(&C, CC_CINTS, 10, 10, 10, 10, 0, "C <ia||jb>", 0, outfile);
  dpd_copy(&C, CC_HBAR, "Wmbej", 0, outfile);
  dpd_buf_close(&C);
  dpd_buf_init(&Wmbej,CC_HBAR, 10, 10, 10, 10, 0, "Wmbej", 0, outfile);
  dpd_scm(&Wmbej, -1, 0, outfile);
  dpd_buf_close(&Wmbej);

  dpd_buf_init(&D, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_copy(&D, CC_TMP0, "WMbEj", 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_init(&WMbEj, CC_TMP0, 0, 5, 0, 5, 0, "WMbEj", 0, outfile);
  dpd_swap24(&WMbEj, CC_HBAR, 10, 11, "WMbEj", 0, outfile);
  dpd_buf_close(&WMbEj);

  dpd_buf_init(&D, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_copy(&D, CC_TMP0, "WmBeJ", 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_init(&WmBeJ, CC_TMP0, 0, 5, 0, 5, 0, "WmBeJ", 0, outfile);
  dpd_swap24(&WmBeJ, CC_HBAR, 10, 11, "WmBeJ", 0, outfile);
  dpd_buf_close(&WmBeJ);

  dpd_buf_init(&C, CC_CINTS, 10, 10, 10, 10, 0, "C <ia|jb>", 0, outfile);
  dpd_copy(&C, CC_HBAR, "WmBEj", 0, outfile);
  dpd_buf_close(&C);
  dpd_buf_init(&WmBEj, CC_HBAR, 10, 10, 10, 10, 0, "WmBEj", 0, outfile);
  dpd_scm(&WmBEj, -1, 0,outfile);
  dpd_buf_close(&WmBEj);

  dpd_buf_init(&C, CC_CINTS, 10, 10, 10, 10, 0, "C <ia|jb>", 0, outfile);
  dpd_copy(&C, CC_HBAR, "WMbeJ", 0, outfile);
  dpd_buf_close(&C);
  dpd_buf_init(&WMbeJ, CC_HBAR, 10, 10, 10, 10, 0, "WMbeJ", 0, outfile);
  dpd_scm(&WMbeJ, -1, 0, outfile);
  dpd_buf_close(&WMbeJ);


  dpd_oe_file_init(&tIA, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_init(&tia, CC_OEI, 0, 1, "tia", 0, outfile);
  
  dpd_buf_init(&F, CC_FINTS, 10, 5, 10, 5, 1, "F <ia|bc>", 0, outfile);
  dpd_buf_init(&WMBEJ, CC_HBAR, 10, 10, 10, 10, 0, "WMBEJ", 0, outfile);
  dpd_contract212(&tIA, &F, &WMBEJ, 1, 2, 1, -1, 1, 0, outfile);
  dpd_buf_close(&WMBEJ);
  dpd_buf_init(&Wmbej, CC_HBAR, 10, 10, 10, 10, 0, "Wmbej", 0, outfile);
  dpd_contract212(&tia, &F, &Wmbej, 1, 2, 1, -1, 1, 0, outfile);
  dpd_buf_close(&Wmbej);
  dpd_buf_close(&F);

  dpd_buf_init(&F, CC_FINTS, 10, 5, 10, 5, 0, "F <ia|bc>", 0, outfile);
  dpd_buf_init(&WMbEj, CC_HBAR, 10, 11, 10, 11, 0, "WMbEj", 0, outfile);
  dpd_contract221(&F, &tia, &WMbEj, 3, 1, 0, 1, 1, 0, outfile);
  dpd_buf_close(&WMbEj);
  dpd_buf_init(&WmBeJ, CC_HBAR, 10, 11, 10, 11, 0, "WmBeJ", 0, outfile);
  dpd_contract221(&F, &tIA, &WmBeJ, 3, 1, 0, 1, 1, 0, outfile);
  dpd_buf_close(&WmBeJ);
  dpd_buf_init(&WMbeJ, CC_HBAR, 10, 10, 10, 10, 0, "WMbeJ", 0, outfile);
  dpd_contract212(&tIA, &F, &WMbeJ, 1, 2, 1, -1, 1, 0, outfile);
  dpd_buf_close(&WMbeJ);
  dpd_buf_init(&WmBEj, CC_HBAR, 10, 10, 10, 10, 0, "WmBEj", 0, outfile);
  dpd_contract212(&tia, &F, &WmBEj, 1, 2, 1, -1, 1, 0, outfile);
  dpd_buf_close(&WmBEj);
  dpd_buf_close(&F);


  dpd_buf_init(&E, CC_EINTS, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)",
	       0, outfile);
  dpd_buf_init(&WMBEJ, CC_HBAR, 10, 10, 10, 10, 0, "WMBEJ", 0, outfile);
  dpd_contract221(&E, &tIA, &WMBEJ, 1, 0, 1, 1, 1, 0, outfile);
  dpd_buf_close(&WMBEJ);
  dpd_buf_init(&Wmbej, CC_HBAR, 10, 10, 10, 10, 0, "Wmbej", 0, outfile);
  dpd_contract221(&E, &tia, &Wmbej, 1, 0, 1, 1, 1, 0, outfile);
  dpd_buf_close(&Wmbej);
  dpd_buf_close(&E);

  dpd_buf_init(&E, CC_EINTS, 11, 0, 11, 0, 0, "E <ai|jk>", 0, outfile);
  dpd_buf_init(&WMbEj, CC_HBAR, 10, 11, 10, 11, 0, "WMbEj", 0, outfile);
  dpd_contract221(&E, &tia, &WMbEj, 3, 0, 1, -1, 1, 0, outfile);
  dpd_buf_close(&WMbEj);
  dpd_buf_init(&WmBeJ, CC_HBAR, 10, 11, 10, 11, 0, "WmBeJ", 0, outfile);
  dpd_contract221(&E, &tIA, &WmBeJ, 3, 0, 1, -1, 1, 0, outfile);
  dpd_buf_close(&WmBeJ);
  dpd_buf_close(&E);

  dpd_buf_init(&E, CC_EINTS, 0, 10, 0, 10, 0, "E <ij|ka>", 0, outfile);
  dpd_buf_init(&WMbeJ, CC_HBAR, 10, 10, 10, 10, 0, "WMbeJ", 0, outfile);
  dpd_contract221(&E, &tia, &WMbeJ, 1, 0, 1, 1, 1, 0, outfile);
  dpd_buf_close(&WMbeJ);
  dpd_buf_init(&WmBEj, CC_HBAR, 10, 10, 10, 10, 0, "WmBEj", 0, outfile);
  dpd_contract221(&E, &tIA, &WmBEj, 1, 0, 1, 1, 1, 0, outfile);
  dpd_buf_close(&WmBEj);
  dpd_buf_close(&E);

  dpd_oe_file_close(&tIA);  dpd_oe_file_close(&tia);


  /* Convert from (MB,JE) to (ME,JB) for remaining terms */

  dpd_buf_init(&WMBEJ, CC_HBAR, 10, 10, 10, 10, 0, "WMBEJ", 0, outfile);
  dpd_swap24(&WMBEJ, CC_TMP0, 10, 10, "WMBEJ", 0, outfile);
  dpd_buf_close(&WMBEJ);
  dpd_buf_init(&WMBEJ, CC_TMP0, 10, 10, 10, 10, 0, "WMBEJ", 0, outfile);
  dpd_copy(&WMBEJ, CC_HBAR, "WMBEJ", 0, outfile);
  dpd_buf_close(&WMBEJ);

  dpd_buf_init(&Wmbej, CC_HBAR, 10, 10, 10, 10, 0, "Wmbej", 0, outfile);
  dpd_swap24(&Wmbej, CC_TMP0, 10, 10, "Wmbej", 0, outfile);
  dpd_buf_close(&Wmbej);
  dpd_buf_init(&Wmbej, CC_TMP0, 10, 10, 10, 10, 0, "Wmbej", 0, outfile);
  dpd_copy(&Wmbej, CC_HBAR, "Wmbej", 0, outfile);
  dpd_buf_close(&Wmbej);

  dpd_buf_init(&WMbEj, CC_HBAR, 10, 11, 10, 11, 0, "WMbEj", 0, outfile);
  dpd_swap34(&WMbEj, CC_TMP0, 10, 10, "WMbEj", 0, outfile);
  dpd_buf_close(&WMbEj);
  dpd_buf_init(&WMbEj, CC_TMP0, 10, 10, 10, 10, 0, "WMbEj", 0, outfile);
  dpd_swap24(&WMbEj, CC_HBAR, 10, 10, "WMbEj", 0, outfile);
  dpd_buf_close(&WMbEj);

  dpd_buf_init(&WmBeJ, CC_HBAR, 10, 11, 10, 11, 0, "WmBeJ", 0, outfile);
  dpd_swap34(&WmBeJ, CC_TMP0, 10, 10, "WmBeJ", 0, outfile);
  dpd_buf_close(&WmBeJ);
  dpd_buf_init(&WmBeJ, CC_TMP0, 10, 10, 10, 10, 0, "WmBeJ", 0, outfile);
  dpd_swap24(&WmBeJ, CC_HBAR, 10, 10, "WmBeJ", 0, outfile);
  dpd_buf_close(&WmBeJ);

  dpd_buf_init(&WMbeJ, CC_HBAR, 10, 10, 10, 10, 0, "WMbeJ", 0, outfile);
  dpd_swap24(&WMbeJ, CC_TMP0, 10, 10, "WMbeJ", 0, outfile);
  dpd_buf_close(&WMbeJ);
  dpd_buf_init(&WMbeJ, CC_TMP0, 10, 10, 10, 10, 0, "WMbeJ", 0, outfile);
  dpd_copy(&WMbeJ, CC_HBAR, "WMbeJ", 0, outfile);
  dpd_buf_close(&WMbeJ);

  dpd_buf_init(&WmBEj, CC_HBAR, 10, 10, 10, 10, 0, "WmBEj", 0, outfile);
  dpd_swap24(&WmBEj, CC_TMP0, 10, 10, "WmBEj", 0, outfile);
  dpd_buf_close(&WmBEj);
  dpd_buf_init(&WmBEj, CC_TMP0, 10, 10, 10, 10, 0, "WmBEj", 0, outfile);
  dpd_copy(&WmBEj, CC_HBAR, "WmBEj", 0, outfile);
  dpd_buf_close(&WmBEj);


  
  dpd_buf_init(&WMBEJ, CC_HBAR, 10, 10, 10, 10, 0, "WMBEJ", 0, outfile);
  dpd_buf_init(&X, CC_MISC, 10, 10, 10, 10, 0, "XIAJB", 0, outfile);
  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)",
	       0, outfile);
  dpd_contract222(&D, &X, &WMBEJ, 0, 0, 1, 1, 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_close(&X);
  dpd_buf_init(&tIAjb, CC_TAMPS, 10, 10, 10, 10, 0, "tIAjb", 0, outfile);
  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)",
	       0, outfile);
  dpd_contract222(&D, &tIAjb, &WMBEJ, 0, 0, 0.5, 1, 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_close(&tIAjb);
  dpd_buf_close(&WMBEJ);
  
  dpd_buf_init(&Wmbej, CC_HBAR, 10, 10, 10, 10, 0, "Wmbej", 0, outfile);
  dpd_buf_init(&X, CC_MISC, 10, 10, 10, 10, 0, "Xiajb", 0, outfile);
  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)",
	       0, outfile);
  dpd_contract222(&D, &X, &Wmbej, 0, 0, 1, 1, 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_close(&X);
  dpd_buf_init(&tiaJB, CC_TAMPS, 10, 10, 10, 10, 0, "tiaJB", 0, outfile);
  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)",
	       0, outfile);
  dpd_contract222(&D, &tiaJB, &Wmbej, 0, 0, 0.5, 1, 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_close(&tiaJB);
  dpd_buf_close(&Wmbej);

  dpd_buf_init(&WMbEj, CC_HBAR, 10, 10, 10, 10, 0, "WMbEj", 0, outfile);
  dpd_buf_init(&tiaJB, CC_TAMPS, 10, 10, 10, 10, 0, "tiaJB", 0, outfile);
  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)",
	       0, outfile);
  dpd_contract222(&D, &tiaJB, &WMbEj, 0, 0, 0.5, 1, 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_close(&tiaJB);
  dpd_buf_init(&X, CC_MISC, 10, 10, 10, 10, 0, "Xiajb", 0, outfile);
  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)",
	       0, outfile);
  dpd_contract222(&D, &X, &WMbEj, 0, 0, 1, 1, 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_close(&X);
  dpd_buf_close(&WMbEj);

  dpd_buf_init(&WmBeJ, CC_HBAR, 10, 10, 10, 10, 0, "WmBeJ", 0, outfile);
  dpd_buf_init(&tIAjb, CC_TAMPS, 10, 10, 10, 10, 0, "tIAjb", 0, outfile);
  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)",
	       0, outfile);
  dpd_contract222(&D, &tIAjb, &WmBeJ, 0, 0, 0.5, 1, 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_close(&tIAjb);
  dpd_buf_init(&X, CC_MISC, 10, 10, 10, 10, 0, "XIAJB", 0, outfile);
  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)", 0, outfile);
  dpd_contract222(&D, &X, &WmBeJ, 0, 0, 1, 1, 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_close(&X);
  dpd_buf_close(&WmBeJ);

  dpd_buf_init(&WMbeJ, CC_HBAR, 10, 10, 10, 10, 0, "WMbeJ", 0, outfile);
  dpd_buf_init(&X, CC_MISC, 10, 10, 10, 10, 0, "XIajB", 0, outfile);
  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)", 0, outfile);
  dpd_contract222(&D, &X, &WMbeJ, 0, 0, 1, 1, 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_close(&X);
  dpd_buf_close(&WMbeJ);

  dpd_buf_init(&WmBEj, CC_HBAR, 10, 10, 10, 10, 0, "WmBEj", 0, outfile);
  dpd_buf_init(&X, CC_MISC, 10, 10, 10, 10, 0, "XiAJb", 0, outfile);
  dpd_buf_init(&D, CC_DINTS, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)",
	       0, outfile);
  dpd_contract222(&D, &X, &WmBEj, 0, 0, 1, 1, 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_close(&X);
  dpd_buf_close(&WmBEj);

/*  timer_off("Wmbej_build", outfile); */
}
