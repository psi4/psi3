#include <stdio.h>
#include <dpd.h>
#include <qt.h>
#define EXTERN
#include "globals.h"

void CT2(void)
{
  struct oe_dpdfile tIA, tia;
  struct dpdbuf Y, C, D, T2new, T2;

/*  timer_on("CT2", outfile); */

  dpd_oe_file_init(&tIA, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_init(&tia, CC_OEI, 0, 1, "tia", 0, outfile);

  /* T(I,E) * C(EM,JB) --> Y(IM,JB) */
  dpd_buf_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, "Y (IM,JB)", 0, outfile);
  dpd_buf_init(&C, CC_CINTS, 11, 10, 11, 10, 0, "C <ia||jb> (bi,ja)",
	       0, outfile);
  dpd_contract212(&tIA, &C, &Y, 1, 0, 0, 1, 0, 0, outfile);  
  dpd_buf_close(&C);

  /* Y(IM,JB) * T(M,A) -> T2(IA,JB) */
  dpd_buf_init(&T2new, CC_TMP1, 10, 10, 10, 10, 0, "T2 (IA,JB)", 0, outfile);
  dpd_contract221(&Y, &tIA, &T2new, 1, 0, 1, 1, 0, 0, outfile);
  dpd_buf_close(&Y);

  /* P(IJ) T2(IA,JB) */
  dpd_swap13(&T2new, CC_TMP2, 10, 10, "T2 (JA,IB)", 0, outfile);
  /* P(AB) T2(IA,JB) */
  dpd_swap24(&T2new, CC_TMP3, 10, 10, "T2 (IB,JA)", 0, outfile);
  dpd_buf_close(&T2new);
  /* P(IJ) P(AB) T2(IA,JB) */
  dpd_buf_init(&T2new, CC_TMP2, 10, 10, 10, 10, 0, "T2 (JA,IB)", 0, outfile);
  dpd_swap24(&T2new, CC_TMP4, 10, 10, "T2 (JB,IA)", 0, outfile);
  dpd_buf_close(&T2new);

  /* T2(IA,JB) - T2(JA,IB) - T2(IB,JA) + T2(JB,IA) --> T2(IA,JB) */
  dpd_buf_init(&T2new, CC_TMP1, 10, 10, 10, 10, 0, "T2 (IA,JB)", 0, outfile);
  dpd_buf_init(&T2, CC_TMP2, 10, 10, 10, 10, 0, "T2 (JA,IB)", 0, outfile);
  dpd_axpy(&T2, &T2new, -1, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_init(&T2, CC_TMP3, 10, 10, 10, 10, 0, "T2 (IB,JA)", 0, outfile);
  dpd_axpy(&T2, &T2new, -1, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_init(&T2, CC_TMP4, 10, 10, 10, 10, 0, "T2 (JB,IA)", 0, outfile);
  dpd_axpy(&T2, &T2new, 1, 0, outfile);
  dpd_buf_close(&T2);

  /* T2(IA,JB) --> T2(IJ,AB) */
  dpd_swap23(&T2new, CC_TMP0, 0, 5, "T2 (IJ,AB)", 0, outfile);
  dpd_buf_close(&T2new);
  dpd_buf_init(&T2new, CC_TMP0, 0, 5, 0, 5, 0, "T2 (IJ,AB)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 0, 5, 2, 7, 0, "New tIJAB",  0, outfile);
  dpd_axpy(&T2new, &T2, 1, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&T2new);


  /* T(i,e) * C(em,jb) --> Y(im,jb) */
  dpd_buf_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, "Y (im,jb)", 0, outfile);
  dpd_buf_init(&C, CC_CINTS,  11, 10, 11, 10, 0, "C <ia||jb> (bi,ja)",
	       0, outfile);
  dpd_contract212(&tia, &C, &Y, 1, 0, 0, 1, 0, 0, outfile);  
  dpd_buf_close(&C);

  /* Y(im,jb) * T(m,a) --> T2(ia,jb) */
  dpd_buf_init(&T2new, CC_TMP1, 10, 10, 10, 10, 0, "T2 (ia,jb)",
	       0, outfile);
  dpd_contract221(&Y, &tia, &T2new, 1, 0, 1, 1, 0, 0, outfile);
  dpd_buf_close(&Y);

  /* P(ij) T2(ja,ib) */
  dpd_swap13(&T2new, CC_TMP2, 10, 10, "T2 (ja,ib)", 0, outfile);
  /* P(ab) T2(ib,ja) */
  dpd_swap24(&T2new, CC_TMP3, 10, 10, "T2 (ib,ja)", 0, outfile);
  dpd_buf_close(&T2new);
  /* P(ij) P(ab) T2(jb,ia) */
  dpd_buf_init(&T2new, CC_TMP2, 10, 10, 10, 10, 0, "T2 (ja,ib)", 0, outfile);
  dpd_swap24(&T2new, CC_TMP4, 10, 10, "T2 (jb,ia)", 0, outfile);
  dpd_buf_close(&T2new);

  /* T2(ia,jb) - T2(ja,ib) - T2(ib,ja) + T2(jb,ia) --> T2(ia,jb) */
  dpd_buf_init(&T2new, CC_TMP1, 10, 10, 10, 10, 0, "T2 (ia,jb)",
	       0, outfile);
  dpd_buf_init(&T2, CC_TMP2, 10, 10, 10, 10, 0, "T2 (ja,ib)", 0, outfile);
  dpd_axpy(&T2, &T2new, -1, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_init(&T2, CC_TMP3, 10, 10, 10, 10, 0, "T2 (ib,ja)", 0, outfile);
  dpd_axpy(&T2, &T2new, -1, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_init(&T2, CC_TMP4, 10, 10, 10, 10, 0, "T2 (jb,ia)", 0, outfile);
  dpd_axpy(&T2, &T2new, 1, 0, outfile);
  dpd_buf_close(&T2);

  /* T2(ia,jb) --> T2(ij,ab) */
  dpd_swap23(&T2new, CC_TMP0, 0, 5, "T2 (ij,ab)", 0, outfile);
  dpd_buf_close(&T2new);

  dpd_buf_init(&T2new, CC_TMP0, 0, 5, 0, 5, 0, "T2 (ij,ab)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 0, 5, 2, 7, 0, "New tijab", 0, outfile);
  dpd_axpy(&T2new, &T2, 1, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&T2new);



  /* T(I,E) * D(EM,jb) --> Y(IM,jb) */
  dpd_buf_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, "Y (IM,jb)", 0, outfile);
  dpd_buf_init(&D, CC_DINTS, 11, 10, 11, 10, 0, "D <ij|ab> (ai,jb)",
	       0, outfile);
  dpd_contract212(&tIA, &D, &Y, 1, 0, 0, -1, 0, 0, outfile);
  dpd_buf_close(&D);

  /* Y(IM,jb) * T(M,A) --> T2(IA,jb) */
  dpd_buf_init(&T2new, CC_TMP1, 10, 10, 10, 10, 0, "T2 (IA,jb)", 0, outfile);
  dpd_contract221(&Y, &tIA, &T2new, 1, 0, 1, 1, 0, 0, outfile);
  dpd_buf_close(&Y);

  /* T2(IA,jb) --> T2(Ij,Ab) */
  dpd_swap23(&T2new, CC_TMP0, 0, 5, "T2 (Ij,Ab)", 0, outfile);
  dpd_buf_close(&T2new);
  dpd_buf_init(&T2new, CC_TMP0, 0, 5, 0, 5, 0, "T2 (Ij,Ab)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 0, 5, 0, 5, 0, "New tIjAb", 0, outfile);
  dpd_axpy(&T2new, &T2, 1, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&T2new);

  /* T(j,e) * D(em,IA) --> Y(jm,IA) */
  dpd_buf_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, "Y (jm,IA)", 0, outfile);
  dpd_buf_init(&D, CC_DINTS, 11, 10, 11, 10, 0, "D <ij|ab> (ai,jb)",
	       0, outfile);
  dpd_contract212(&tia, &D, &Y, 1, 0, 0, -1, 0, 0, outfile);
  dpd_buf_close(&D);

  /* Y(jm,IA) * T(m,b) --> T2(jb,IA) */
  dpd_buf_init(&T2new, CC_TMP1, 10, 10, 10, 10, 0, "T2 (jb,IA)", 0, outfile);
  dpd_contract221(&Y, &tia, &T2new, 1, 0, 1, 1, 0, 0, outfile);
  dpd_buf_close(&Y);
  
  /* T2(jb,IA) --> T2(jI,bA) */
  dpd_swap23(&T2new, CC_TMP0, 0, 5, "T2 (jI,bA)", 0, outfile);
  dpd_buf_close(&T2new);
  dpd_buf_init(&T2new, CC_TMP0, 0, 5, 0, 5, 0, "T2 (jI,bA)", 0, outfile);
  /* T2(jI,bA) --> T2(Ij,bA) */
  dpd_swap12(&T2new, CC_TMP1, 0, 5, "T2 (Ij,bA)", 0, outfile);
  dpd_buf_close(&T2new);
  /* T2(Ij,bA) --> T2(Ij,Ab) */
  dpd_buf_init(&T2new, CC_TMP1, 0, 5, 0, 5, 0, "T2 (Ij,bA)", 0, outfile);
  dpd_swap34(&T2new, CC_TMP0, 0, 5, "T2 (Ij,Ab)", 0, outfile);
  dpd_buf_close(&T2new);
  dpd_buf_init(&T2new, CC_TMP0, 0, 5, 0, 5, 0, "T2 (Ij,Ab)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 0, 5, 0, 5, 0, "New tIjAb", 0, outfile);
  dpd_axpy(&T2new, &T2, 1, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&T2new);



  /* T(I,E) * C(Em,jA) --> Y(Im,jA) */
  dpd_buf_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, "Y (Im,jA)", 0, outfile);
  dpd_buf_init(&C, CC_CINTS, 11, 10, 11, 10, 0, "C <ia|jb> (bi,ja)", 0, outfile);
  dpd_contract212(&tIA, &C, &Y, 1, 0, 0, -1, 0, 0, outfile);
  dpd_buf_close(&C);

  /* Y(Im,jA) * T(m,b) --> T2(Ib,jA) */
  dpd_buf_init(&T2new, CC_TMP1, 10, 10, 10, 10, 0, "T2 (Ib,jA)", 0, outfile);
  dpd_contract221(&Y, &tia, &T2new, 1, 0, 1, 1, 0, 0, outfile);
  dpd_buf_close(&Y);
  /* T2(Ib,jA) --> T2(Ij,bA) */
  dpd_swap23(&T2new, CC_TMP0, 0, 5, "T2 (Ij,bA)", 0, outfile);
  dpd_buf_close(&T2new);
  /* T2(Ij,bA) --> T2(Ij,Ab) */
  dpd_buf_init(&T2new, CC_TMP0, 0, 5, 0, 5, 0, "T2 (Ij,bA)", 0, outfile);
  dpd_swap34(&T2new, CC_TMP1, 0, 5, "T2 (Ij,Ab)", 0, outfile);
  dpd_buf_close(&T2new);
  dpd_buf_init(&T2new, CC_TMP1, 0, 5, 0, 5, 0, "T2 (Ij,Ab)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 0, 5, 0, 5, 0, "New tIjAb", 0, outfile);
  dpd_axpy(&T2new, &T2, 1, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&T2new);


  /* T(j,e) * C(eM,Ib) --> Y(jM,Ib) */
  dpd_buf_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, "Y (jM,Ib)", 0, outfile);
  dpd_buf_init(&C, CC_CINTS, 11, 10, 11, 10, 0, "C <ia|jb> (bi,ja)", 0, outfile);
  dpd_contract212(&tia, &C, &Y, 1, 0, 0, -1, 0, 0, outfile);
  dpd_buf_close(&C);

  /* Y(jM,Ib) * T(M,A) --> T2(jA,Ib) */
  dpd_buf_init(&T2new, CC_TMP1, 10, 10, 10, 10, 0, "T2 (jA,Ib)", 0, outfile);
  dpd_contract221(&Y, &tIA, &T2new, 1, 0, 1, 1, 0, 0, outfile);
  dpd_buf_close(&Y);
  /* T2(jA,Ib) --> T2(jI,Ab) */
  dpd_swap23(&T2new, CC_TMP0, 0, 5, "T2 (jI,Ab)", 0, outfile);
  dpd_buf_close(&T2new);
  /* T2(jI,Ab) --> T2(Ij,Ab) */
  dpd_buf_init(&T2new, CC_TMP0, 0, 5, 0, 5, 0, "T2 (jI,Ab)", 0, outfile);
  dpd_swap12(&T2new, CC_TMP1, 0, 5, "T2 (Ij,Ab)", 0, outfile);
  dpd_buf_close(&T2new);
  dpd_buf_init(&T2new, CC_TMP1, 0, 5, 0, 5, 0, "T2 (Ij,Ab)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 0, 5, 0, 5, 0, "New tIjAb", 0, outfile);
  dpd_axpy(&T2new, &T2, 1, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&T2new);

  dpd_oe_file_close(&tIA); dpd_oe_file_close(&tia);

/*  timer_off("CT2", outfile); */
}
