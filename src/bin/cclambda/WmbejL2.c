#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void WmbejL2(void)
{
  struct dpdbuf newL2, L2, W, Z, Z2;

  /* RHS += P(ij)P(ab)Limae * Wjebm */

  dpd_buf_init(&Z, CC_TMP0, 10, 10, 10, 10, 0, "Z(IA,JB)", 0, outfile);
  dpd_buf_init(&L2, CC_LAMPS, 10, 10, 10, 10, 0, "LIAJB", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 10, 10, 10, 0, "WMBEJ", 0, outfile);
  dpd_contract222(&L2, &W, &Z, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&L2);
  dpd_buf_init(&L2, CC_LAMPS, 10, 10, 10, 10, 0, "LIAjb", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 10, 10, 10, 0, "WMbEj", 0, outfile);
  dpd_contract222(&L2, &W, &Z, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&L2);
  dpd_swap13(&Z, CC_TMP1, 10, 10, "Z(JA,IB)", 0, outfile);
  dpd_swap24(&Z, CC_TMP2, 10, 10, "Z(IB,JA)", 0, outfile);
  dpd_swapbk(&Z, CC_TMP3, 10, 10, "Z(JB,IA)", 0, outfile);
  dpd_buf_init(&Z2, CC_TMP1, 10, 10, 10, 10, 0, "Z(JA,IB)", 0, outfile);
  dpd_axpy(&Z2, &Z, -1.0, 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_init(&Z2, CC_TMP2, 10, 10, 10, 10, 0, "Z(IB,JA)", 0, outfile);
  dpd_axpy(&Z2, &Z, -1.0, 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_init(&Z2, CC_TMP3, 10, 10, 10, 10, 0, "Z(JB,IA)", 0, outfile);
  dpd_axpy(&Z2, &Z, 1.0, 0, outfile);
  dpd_buf_close(&Z2);
  dpd_swap23(&Z, CC_TMP1, 0, 5, "Z(IJ,AB)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&Z, CC_TMP1, 0, 5, 0, 5, 0, "Z(IJ,AB)", 0, outfile);
  dpd_buf_init(&newL2, CC_LAMPS, 0, 5, 2, 7, 0, "New LIJAB", 0, outfile);
  dpd_axpy(&Z, &newL2, 1.0, 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&newL2);

  dpd_buf_init(&Z, CC_TMP0, 10, 10, 10, 10, 0, "Z(ia,jb)", 0, outfile);
  dpd_buf_init(&L2, CC_LAMPS, 10, 10, 10, 10, 0, "Liajb", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 10, 10, 10, 0, "Wmbej", 0, outfile);
  dpd_contract222(&L2, &W, &Z, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&L2);
  dpd_buf_init(&L2, CC_LAMPS, 10, 10, 10, 10, 0, "LiaJB", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 10, 10, 10, 0, "WmBeJ", 0, outfile);
  dpd_contract222(&L2, &W, &Z, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&L2);
  dpd_swap13(&Z, CC_TMP1, 10, 10, "Z(ja,ib)", 0, outfile);
  dpd_swap24(&Z, CC_TMP2, 10, 10, "Z(ib,ja)", 0, outfile);
  dpd_swapbk(&Z, CC_TMP3, 10, 10, "Z(jb,ia)", 0, outfile);
  dpd_buf_init(&Z2, CC_TMP1, 10, 10, 10, 10, 0, "Z(ja,ib)", 0, outfile);
  dpd_axpy(&Z2, &Z, -1.0, 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_init(&Z2, CC_TMP2, 10, 10, 10, 10, 0, "Z(ib,ja)", 0, outfile);
  dpd_axpy(&Z2, &Z, -1.0, 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_init(&Z2, CC_TMP3, 10, 10, 10, 10, 0, "Z(jb,ia)", 0, outfile);
  dpd_axpy(&Z2, &Z, 1.0, 0, outfile);
  dpd_buf_close(&Z2);
  dpd_swap23(&Z, CC_TMP1, 0, 5, "Z(ij,ab)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&Z, CC_TMP1, 0, 5, 0, 5, 0, "Z(ij,ab)", 0, outfile);
  dpd_buf_init(&newL2, CC_LAMPS, 0, 5, 2, 7, 0, "New Lijab", 0, outfile);
  dpd_axpy(&Z, &newL2, 1.0, 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&newL2);


  dpd_buf_init(&Z, CC_TMP0, 10, 10, 10, 10, 0, "Z(IA,jb)", 0, outfile);
  dpd_buf_init(&L2, CC_LAMPS, 10, 10, 10, 10, 0, "LIAJB", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 10, 10, 10, 0, "WmBeJ", 0, outfile);
  dpd_contract222(&L2, &W, &Z, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&L2);
  dpd_buf_close(&W);
  dpd_buf_init(&L2, CC_LAMPS, 10, 10, 10, 10, 0, "LIAjb", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 10, 10, 10, 0, "Wmbej", 0, outfile);
  dpd_contract222(&L2, &W, &Z, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&L2);
  dpd_buf_init(&L2, CC_LAMPS, 10, 10, 10, 10, 0, "Liajb", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 10, 10, 10, 0, "WMbEj", 0, outfile);
  dpd_contract222(&W, &L2, &Z, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&L2);
  dpd_buf_close(&W);
  dpd_buf_init(&L2, CC_LAMPS, 10, 10, 10, 10, 0, "LiaJB", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 10, 10, 10, 0, "WMBEJ", 0, outfile);
  dpd_contract222(&W, &L2, &Z, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&L2);
  dpd_swap23(&Z, CC_TMP1, 0, 5, "Z(Ij,Ab)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&Z, CC_TMP1, 0, 5, 0, 5, 0, "Z(Ij,Ab)", 0, outfile);
  dpd_buf_init(&newL2, CC_LAMPS, 0, 5, 0, 5, 0, "New LIjAb", 0, outfile);
  dpd_axpy(&Z, &newL2, 1.0, 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&newL2);

  dpd_buf_init(&Z, CC_TMP0, 10, 10, 10, 10, 0, "Z(Ib,jA)", 0, outfile);
  dpd_buf_init(&L2, CC_LAMPS, 10, 10, 10, 10, 0, "LIbjA", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 10, 10, 10, 0, "WMbeJ", 0, outfile);
  dpd_contract222(&W, &L2, &Z, 0, 1, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&L2);
  dpd_buf_close(&W);
  dpd_buf_init(&L2, CC_LAMPS, 10, 10, 10, 10, 0, "LjAIb", 0, outfile);
  dpd_buf_init(&W, CC_HBAR, 10, 10, 10, 10, 0, "WmBEj", 0, outfile);
  dpd_contract222(&L2, &W, &Z, 1, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&L2);
  dpd_buf_close(&W);
  dpd_swap23(&Z, CC_TMP1, 0, 5, "Z(Ij,bA)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&Z, CC_TMP1, 0, 5, 0, 5, 0, "Z(Ij,bA)", 0, outfile);
  dpd_swap34(&Z, CC_TMP0, 0, 5, "Z(Ij,Ab)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&Z, CC_TMP0, 0, 5, 0, 5, 0, "Z(Ij,Ab)", 0, outfile);
  dpd_buf_init(&newL2, CC_LAMPS, 0, 5, 0, 5, 0, "New LIjAb", 0, outfile);
  dpd_axpy(&Z, &newL2, 1.0, 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&newL2);
}
