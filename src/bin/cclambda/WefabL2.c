#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void WefabL2(void)
{
  struct dpdbuf Lijab, LIJAB, LIjAb;
  struct dpdbuf newLijab, newLIJAB, newLIjAb;
  struct dpdbuf Tau, Z, L2, B, D, F, Ltmp;
  struct oe_dpdfile tIA, tia;

  /* RHS += Wefab*Lijef  */
  dpd_oe_file_init(&tIA, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_init(&tia, CC_OEI, 0, 1, "tia", 0, outfile);

  dpd_buf_init(&newLIJAB, CC_LAMPS, 2, 7, 2, 7, 0, "New LIJAB",
              0, outfile);

  dpd_buf_init(&LIJAB, CC_LAMPS, 2, 7, 2, 7, 0, "LIJAB", 0, outfile);
  dpd_buf_init(&B, CC_BINTS, 7, 7, 5, 5, 1, "B <ab|cd>", 0, outfile);
  dpd_contract222(&LIJAB, &B, &newLIJAB, 0, 1, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&B);
  dpd_buf_close(&LIJAB);

  dpd_buf_init(&LIJAB, CC_LAMPS, 2, 5, 2, 7, 0, "LIJAB", 0, outfile);
  dpd_buf_init(&Ltmp, CC_TMP0, 2, 10, 2, 10, 0, "Ltmp (I>J,MF)", 0, outfile);
  dpd_buf_init(&F, CC_FINTS, 10, 7, 10, 5, 1, "F <ia|bc>", 0, outfile);
  dpd_contract212(&tIA, &LIJAB, &Ltmp, 1, 2, 1, 1.0, 0.0, 0, outfile);
  dpd_contract222(&Ltmp, &F, &newLIJAB, 0, 1, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&F);
  dpd_buf_close(&Ltmp);
  dpd_buf_close(&LIJAB);

  dpd_buf_init(&Z, CC_TMP0, 2, 2, 2, 2, 0, "Z(IJ,MN)", 0, outfile);
  dpd_buf_init(&Tau, CC_TAMPS, 2, 7, 2, 7, 0, "tauIJAB", 0, outfile);
  dpd_buf_init(&L2, CC_LAMPS, 2, 7, 2, 7, 0, "LIJAB", 0, outfile);
  dpd_contract222(&L2, &Tau, &Z, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&L2);
  dpd_buf_close(&Tau);
  dpd_buf_init(&D, CC_DINTS, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)",
	       0, outfile);
  dpd_contract222(&Z, &D, &newLIJAB, 0, 1, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_close(&Z);
  dpd_buf_close(&newLIJAB);

  dpd_buf_init(&newLijab, CC_LAMPS, 2, 7, 2, 7, 0, "New Lijab",
              0, outfile);

  dpd_buf_init(&Lijab, CC_LAMPS, 2, 7, 2, 7, 0, "Lijab", 0, outfile);
  dpd_buf_init(&B, CC_BINTS, 7, 7, 5, 5, 1, "B <ab|cd>", 0, outfile);
  dpd_contract222(&Lijab, &B, &newLijab, 0, 1, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&B);
  dpd_buf_close(&Lijab);

  dpd_buf_init(&Lijab, CC_LAMPS, 2, 5, 2, 7, 0, "Lijab", 0, outfile);
  dpd_buf_init(&Ltmp, CC_TMP0, 2, 10, 2, 10, 0, "Ltmp (i>j,mf)", 0, outfile);
  dpd_buf_init(&F, CC_FINTS, 10, 7, 10, 5, 1, "F <ia|bc>", 0, outfile);
  dpd_contract212(&tia, &Lijab, &Ltmp, 1, 2, 1, 1.0, 0.0, 0, outfile);
  dpd_contract222(&Ltmp, &F, &newLijab, 0, 1, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&F);
  dpd_buf_close(&Ltmp);
  dpd_buf_close(&Lijab);

  dpd_buf_init(&Z, CC_TMP0, 2, 2, 2, 2, 0, "Z(ij,mn)", 0, outfile);
  dpd_buf_init(&Tau, CC_TAMPS, 2, 7, 2, 7, 0, "tauijab", 0, outfile);
  dpd_buf_init(&L2, CC_LAMPS, 2, 7, 2, 7, 0, "Lijab", 0, outfile);
  dpd_contract222(&L2, &Tau, &Z, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&L2);
  dpd_buf_close(&Tau);
  dpd_buf_init(&D, CC_DINTS, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)",
	       0, outfile);
  dpd_contract222(&Z, &D, &newLijab, 0, 1, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_close(&Z);
  dpd_buf_close(&newLijab);


  dpd_buf_init(&newLIjAb, CC_LAMPS, 0, 5, 0, 5, 0, "New LIjAb",
             0, outfile);

  dpd_buf_init(&LIjAb, CC_LAMPS, 0, 5, 0, 5, 0, "LIjAb", 0, outfile);

  dpd_buf_init(&B, CC_BINTS, 5, 5, 5, 5, 0, "B <ab|cd>", 0, outfile);
  dpd_contract222(&LIjAb, &B, &newLIjAb, 0, 1, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&B);
  dpd_buf_init(&Ltmp, CC_TMP1, 0, 11, 0, 11, 0, "Lt (Ij,Em)", 0, outfile);
  dpd_contract221(&LIjAb, &tia, &Ltmp, 3, 1, 0, 1.0, 0.0, 0, outfile);
  dpd_swap34(&Ltmp, CC_TMP2, 0, 10, "Lt (Ij,mE)", 0, outfile);
  dpd_buf_close(&Ltmp);
  dpd_buf_init(&Ltmp, CC_TMP3, 0, 10, 0, 10, 0, "Lt (Ij,Mf)", 0, outfile);
  dpd_contract212(&tIA, &LIjAb, &Ltmp, 1, 2, 1, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&Ltmp);

  dpd_buf_close(&LIjAb);

  dpd_buf_init(&F, CC_FINTS, 10, 5, 10, 5, 0, "F <ia|bc>", 0, outfile);
  dpd_buf_init(&Ltmp, CC_TMP3, 0, 10, 0, 10, 0, "Lt (Ij,Mf)", 0, outfile);
  dpd_contract222(&Ltmp, &F, &newLIjAb, 0, 1, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&Ltmp);
  dpd_swap34(&F, CC_TMP0, 10, 5, "<me|ab> (mE,Ab)", 0, outfile);
  dpd_buf_close(&F);

  dpd_buf_init(&F, CC_TMP0, 10, 5, 10, 5, 0, "<me|ab> (mE,Ab)", 0, outfile);
  dpd_buf_init(&Ltmp, CC_TMP2, 0, 10, 0, 10, 0, "Lt (Ij,mE)", 0, outfile);
  dpd_contract222(&Ltmp, &F, &newLIjAb, 0, 1, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&Ltmp);
  dpd_buf_close(&F);

  dpd_buf_init(&Z, CC_TMP0, 0, 0, 0, 0, 0, "Z(Ij,Mn)", 0, outfile);
  dpd_buf_init(&Tau, CC_TAMPS, 0, 5, 0, 5, 0, "tauIjAb", 0, outfile);
  dpd_buf_init(&L2, CC_LAMPS, 0, 5, 0, 5, 0, "LIjAb", 0, outfile);
  dpd_contract222(&L2, &Tau, &Z, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&L2);
  dpd_buf_close(&Tau);
  dpd_buf_init(&D, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_contract222(&Z, &D, &newLIjAb, 0, 1, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&D);
  dpd_buf_close(&Z);

  dpd_buf_close(&newLIjAb);

  dpd_oe_file_close(&tIA);
  dpd_oe_file_close(&tia);

}
