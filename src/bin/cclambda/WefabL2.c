#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void WefabL2(void)
{
  dpdbuf4 Lijab, LIJAB, LIjAb;
  dpdbuf4 newLijab, newLIJAB, newLIjAb;
  dpdbuf4 Tau, Z, L2, B, D, F, Ltmp;
  dpdfile2 tIA, tia;

  /* RHS += Wefab*Lijef  */
  dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

  dpd_buf4_init(&newLIJAB, CC_LAMPS, 0, 2, 7, 2, 7, 0, "New LIJAB");

  if(!params.aobasis) {
    dpd_buf4_init(&LIJAB, CC_LAMPS, 0, 2, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&B, CC_BINTS, 0, 7, 7, 5, 5, 1, "B <ab|cd>");
    dpd_contract444(&LIJAB, &B, &newLIJAB, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&B);
    dpd_buf4_close(&LIJAB);
  }

  dpd_buf4_init(&LIJAB, CC_LAMPS, 0, 2, 5, 2, 7, 0, "LIJAB");
  dpd_buf4_init(&Ltmp, CC_TMP0, 0, 2, 10, 2, 10, 0, "Ltmp (I>J,MF)");
  dpd_buf4_init(&F, CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
  dpd_contract244(&tIA, &LIJAB, &Ltmp, 1, 2, 1, 1.0, 0.0);
  dpd_contract444(&Ltmp, &F, &newLIJAB, 0, 1, -1.0, 1.0);
  dpd_buf4_close(&F);
  dpd_buf4_close(&Ltmp);
  dpd_buf4_close(&LIJAB);

  dpd_buf4_init(&Z, CC_TMP0, 0, 2, 2, 2, 2, 0, "Z(IJ,MN)");
  dpd_buf4_init(&Tau, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
  dpd_buf4_init(&L2, CC_LAMPS, 0, 2, 7, 2, 7, 0, "LIJAB");
  dpd_contract444(&L2, &Tau, &Z, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&L2);
  dpd_buf4_close(&Tau);
  dpd_buf4_init(&D, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
  dpd_contract444(&Z, &D, &newLIJAB, 0, 1, 1.0, 1.0);
  dpd_buf4_close(&D);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&newLIJAB);

  dpd_buf4_init(&newLijab, CC_LAMPS, 0, 2, 7, 2, 7, 0, "New Lijab");

  if(!params.aobasis) {
    dpd_buf4_init(&Lijab, CC_LAMPS, 0, 2, 7, 2, 7, 0, "Lijab");
    dpd_buf4_init(&B, CC_BINTS, 0, 7, 7, 5, 5, 1, "B <ab|cd>");
    dpd_contract444(&Lijab, &B, &newLijab, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&B);
    dpd_buf4_close(&Lijab);
  }

  dpd_buf4_init(&Lijab, CC_LAMPS, 0, 2, 5, 2, 7, 0, "Lijab");
  dpd_buf4_init(&Ltmp, CC_TMP0, 0, 2, 10, 2, 10, 0, "Ltmp (i>j,mf)");
  dpd_buf4_init(&F, CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
  dpd_contract244(&tia, &Lijab, &Ltmp, 1, 2, 1, 1.0, 0.0);
  dpd_contract444(&Ltmp, &F, &newLijab, 0, 1, -1.0, 1.0);
  dpd_buf4_close(&F);
  dpd_buf4_close(&Ltmp);
  dpd_buf4_close(&Lijab);

  dpd_buf4_init(&Z, CC_TMP0, 0, 2, 2, 2, 2, 0, "Z(ij,mn)");
  dpd_buf4_init(&Tau, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
  dpd_buf4_init(&L2, CC_LAMPS, 0, 2, 7, 2, 7, 0, "Lijab");
  dpd_contract444(&L2, &Tau, &Z, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&L2);
  dpd_buf4_close(&Tau);
  dpd_buf4_init(&D, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
  dpd_contract444(&Z, &D, &newLijab, 0, 1, 1.0, 1.0);
  dpd_buf4_close(&D);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&newLijab);


  dpd_buf4_init(&newLIjAb, CC_LAMPS, 0, 0, 5, 0, 5, 0, "New LIjAb");

  dpd_buf4_init(&LIjAb, CC_LAMPS, 0, 0, 5, 0, 5, 0, "LIjAb");

  if(!params.aobasis) {
    dpd_buf4_init(&B, CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
    dpd_contract444(&LIjAb, &B, &newLIjAb, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&B);
  }

  dpd_buf4_init(&Ltmp, CC_TMP1, 0, 0, 11, 0, 11, 0, "Lt (Ij,Em)");
  dpd_contract424(&LIjAb, &tia, &Ltmp, 3, 1, 0, 1.0, 0.0);
  dpd_buf4_sort(&Ltmp, CC_TMP2, pqsr, 0, 10, "Lt (Ij,mE)");
  dpd_buf4_close(&Ltmp);
  dpd_buf4_init(&Ltmp, CC_TMP3, 0, 0, 10, 0, 10, 0, "Lt (Ij,Mf)");
  dpd_contract244(&tIA, &LIjAb, &Ltmp, 1, 2, 1, 1.0, 0.0);
  dpd_buf4_close(&Ltmp);

  dpd_buf4_close(&LIjAb);

  dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_buf4_init(&Ltmp, CC_TMP3, 0, 0, 10, 0, 10, 0, "Lt (Ij,Mf)");
  dpd_contract444(&Ltmp, &F, &newLIjAb, 0, 1, -1.0, 1.0);
  dpd_buf4_close(&Ltmp);
  dpd_buf4_sort(&F, CC_TMP0, pqsr, 10, 5, "<me|ab> (mE,Ab)");
  dpd_buf4_close(&F);

  dpd_buf4_init(&F, CC_TMP0, 0, 10, 5, 10, 5, 0, "<me|ab> (mE,Ab)");
  dpd_buf4_init(&Ltmp, CC_TMP2, 0, 0, 10, 0, 10, 0, "Lt (Ij,mE)");
  dpd_contract444(&Ltmp, &F, &newLIjAb, 0, 1, -1.0, 1.0);
  dpd_buf4_close(&Ltmp);
  dpd_buf4_close(&F);

  dpd_buf4_init(&Z, CC_TMP0, 0, 0, 0, 0, 0, 0, "Z(Ij,Mn)");
  dpd_buf4_init(&Tau, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
  dpd_buf4_init(&L2, CC_LAMPS, 0, 0, 5, 0, 5, 0, "LIjAb");
  dpd_contract444(&L2, &Tau, &Z, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&L2);
  dpd_buf4_close(&Tau);
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_contract444(&Z, &D, &newLIjAb, 0, 1, 1.0, 1.0);
  dpd_buf4_close(&D);
  dpd_buf4_close(&Z);

  dpd_buf4_close(&newLIjAb);

  dpd_file2_close(&tIA);
  dpd_file2_close(&tia);

}
