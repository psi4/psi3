#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/** The RHF/ROHF version need a significant amount of work.  The 
original version of this code invoved some strange things with
contractions of amplitudes. **/

void WefabL2(void)
{
  dpdbuf4 Lijab, LIJAB, LIjAb;
  dpdbuf4 newLijab, newLIJAB, newLIjAb;
  dpdbuf4 Tau, T2, Z, Z1, Z2, L, L2, B, D, F, Ltmp;
  dpdfile2 tIA, tia;

  /* RHS += Wefab*Lijef  */
  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    dpd_buf4_init(&newLIJAB, CC_LAMPS, 0, 2, 7, 2, 7, 0, "New LIJAB");

    dpd_buf4_init(&LIJAB, CC_LAMPS, 0, 2, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&B, CC_BINTS, 0, 7, 7, 5, 5, 1, "B <ab|cd>");
    dpd_contract444(&LIJAB, &B, &newLIJAB, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&B);
    dpd_buf4_close(&LIJAB);

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

    dpd_buf4_init(&Lijab, CC_LAMPS, 0, 2, 7, 2, 7, 0, "Lijab");
    dpd_buf4_init(&B, CC_BINTS, 0, 7, 7, 5, 5, 1, "B <ab|cd>");
    dpd_contract444(&Lijab, &B, &newLijab, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&B);
    dpd_buf4_close(&Lijab);

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

    dpd_buf4_init(&B, CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
    dpd_contract444(&LIjAb, &B, &newLIjAb, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&B);

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
  else if(params.ref == 2) { /** UHF **/

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");


    /** Z(AB,IJ) = <AB||CD> L(IJ,CD) **/
    dpd_buf4_init(&Z, CC_TMP1, 0, 7, 2, 7, 2, 0, "Z(AB,IJ)");
    dpd_buf4_init(&LIJAB, CC_LAMPS, 0, 2, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&B, CC_BINTS, 0, 7, 7, 5, 5, 1, "B <AB|CD>");
    dpd_contract444(&B, &LIJAB, &Z, 0, 0, 1, 0);
    dpd_buf4_close(&B);
    dpd_buf4_close(&LIJAB);
    /** Z(AB,IJ) --> New L(IJ,AB) **/
    dpd_buf4_sort_axpy(&Z, CC_LAMPS, rspq, 2, 7, "New LIJAB", 1);
    dpd_buf4_close(&Z);


    /** Z(MB,IJ) = - <MB||EF> L(IJ,EF) **/
    dpd_buf4_init(&Z, CC_TMP1, 0, 20, 2, 20, 2, 0, "Z(MB,IJ)");
    dpd_buf4_init(&LIJAB, CC_LAMPS, 0, 2, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&F, CC_FINTS, 0, 20, 7, 20, 5, 1, "F <IA|BC>");
    dpd_contract444(&F, &LIJAB, &Z, 0, 0, -1, 0);
    dpd_buf4_close(&F);
    dpd_buf4_close(&LIJAB);
    /** Z(IJ,AB) <-- T(M,A) Z(MB,IJ) **/
    dpd_buf4_init(&Z1, CC_TMP1, 0, 2, 5, 2, 5, 0, "Z(IJ,AB)");
    dpd_contract244(&tIA, &Z, &Z1, 0, 0, 1, 1, 0);
    dpd_buf4_close(&Z);
    /** Z(IJ,AB) --> Z'(IJ,BA) **/
    dpd_buf4_sort(&Z1, CC_TMP1, qprs, 2, 5, "Z'(IJ,BA)");
    dpd_buf4_init(&Z2, CC_TMP1, 0, 2, 5, 2, 5, 0, "Z'(IJ,BA)");
    /** Z(IJ,AB) = Z(IJ,AB) - Z'(IJ,BA) **/
    dpd_buf4_axpy(&Z2, &Z1, -1);
    dpd_buf4_close(&Z2);
    /** Z(IJ,AB) --> New L(IJ,AB) **/
    dpd_buf4_init(&newLIJAB, CC_LAMPS, 0, 2, 5, 2, 7, 0, "New LIJAB");
    dpd_buf4_axpy(&Z1, &newLIJAB, 1);
    dpd_buf4_close(&newLIJAB);
    dpd_buf4_close(&Z1);

    /** Z(IJ,MN) = L(IJ,EF) <MN||EF> **/
    dpd_buf4_init(&Z, CC_TMP1, 0, 2, 2, 2, 2, 0, "Z(IJ,MN)");
    dpd_buf4_init(&D, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
    dpd_buf4_init(&L2, CC_LAMPS, 0, 2, 7, 2, 7, 0, "LIJAB");
    dpd_contract444(&L2, &D, &Z, 0, 0, 1, 0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&D);
    /** New L(IJ,AB) <-- Z(IJ,MN) Tau(MN,AB) **/
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    dpd_buf4_init(&newLIJAB, CC_LAMPS, 0, 2, 7, 2, 7, 0, "New LIJAB");
    dpd_contract444(&Z, &T2, &newLIJAB, 0, 1, 1, 1);
    dpd_buf4_close(&newLIJAB);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&Z);


    /** Z(ab,ij) = <ab||cd> L(ij,cd) **/
    dpd_buf4_init(&Z, CC_TMP1, 0, 17, 12, 17, 12, 0, "Z(ab,ij)");
    dpd_buf4_init(&LIJAB, CC_LAMPS, 0, 12, 17, 12, 17, 0, "Lijab");
    dpd_buf4_init(&B, CC_BINTS, 0, 17, 17, 15, 15, 1, "B <ab|cd>");
    dpd_contract444(&B, &LIJAB, &Z, 0, 0, 1, 0);
    dpd_buf4_close(&B);
    dpd_buf4_close(&LIJAB);
    /** Z(ab,ij) --> New L(ij,ab) **/
    dpd_buf4_sort_axpy(&Z, CC_LAMPS, rspq, 12, 17, "New Lijab", 1);
    dpd_buf4_close(&Z);

    /** Z(mb,ij) = - <mb||ef> L(ij,ef) **/
    dpd_buf4_init(&Z, CC_TMP1, 0, 30, 12, 30, 12, 0, "Z(mb,ij)");
    dpd_buf4_init(&LIJAB, CC_LAMPS, 0, 12, 17, 12, 17, 0, "Lijab");
    dpd_buf4_init(&F, CC_FINTS, 0, 30, 17, 30, 15, 1, "F <ia|bc>");
    dpd_contract444(&F, &LIJAB, &Z, 0, 0, -1, 0);
    dpd_buf4_close(&F);
    dpd_buf4_close(&LIJAB);
    /** Z(ij,ab) <-- T(m,a) Z(mb,ij) **/
    dpd_buf4_init(&Z1, CC_TMP1, 0, 12, 15, 12, 15, 0, "Z(ij,ab)");
    dpd_contract244(&tia, &Z, &Z1, 0, 0, 1, 1, 0);
    dpd_buf4_close(&Z);
    /** Z(ij,ab) --> Z'(ij,ba) **/
    dpd_buf4_sort(&Z1, CC_TMP1, pqsr, 12, 15, "Z'(ij,ba)");
    dpd_buf4_init(&Z2, CC_TMP1, 0, 12, 15, 12, 15, 0, "Z'(ij,ba)");
    /** Z(ij,ab) = Z(ij,ab) - Z'(ij,ba) **/
    dpd_buf4_axpy(&Z2, &Z1, -1);
    dpd_buf4_close(&Z2);
    /** Z(ij,ab) --> New L(ij,ab) **/
    dpd_buf4_init(&newLijab, CC_LAMPS, 0, 12, 15, 12, 17, 0, "New Lijab");
    dpd_buf4_axpy(&Z1, &newLijab, 1);
    dpd_buf4_close(&newLijab);
    dpd_buf4_close(&Z1);

    /** Z(ij,mn) = L(ij,ef) <mn||ef> **/
    dpd_buf4_init(&Z, CC_TMP1, 0, 12, 12, 12, 12, 0, "Z(ij,mn)");
    dpd_buf4_init(&D, CC_DINTS, 0, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
    dpd_buf4_init(&L2, CC_LAMPS, 0, 12, 17, 12, 17, 0, "Lijab");
    dpd_contract444(&L2, &D, &Z, 0, 0, 1, 0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&D);
    /** New L(ij,ab) <-- Z(ij,mn) Tau(mn,ab) **/
    dpd_buf4_init(&newLijab, CC_LAMPS, 0, 12, 17, 12, 17, 0, "New Lijab");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
    dpd_contract444(&Z, &T2, &newLijab, 0, 1, 1, 1);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&newLijab);
    dpd_buf4_close(&Z);





    /** Z(Ab,Ij) = <Ab|Cd> L(Ij,Cd) **/
    dpd_buf4_init(&Z, CC_TMP1, 0, 28, 22, 28, 22, 0, "Z(Ab,Ij)");
    dpd_buf4_init(&LIjAb, CC_LAMPS, 0, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_init(&B, CC_BINTS, 0, 28, 28, 28, 28, 0, "B <Ab|Cd>");
    dpd_contract444(&B, &LIjAb, &Z, 0, 0, 1, 0);
    dpd_buf4_close(&B);
    dpd_buf4_close(&LIjAb);
    /** Z(Ab,Ij) --> New L(Ij,Ab) **/
    dpd_buf4_sort_axpy(&Z, CC_LAMPS, rspq, 22, 28, "New LIjAb", 1);
    dpd_buf4_close(&Z);


    /** Z(Ij,Am) = L(Ij,Fe) F(Fe,Am) **/
    dpd_buf4_init(&Z, CC_TMP1, 0, 22, 26, 22, 26, 0, "Z(Ij,Am)");
    dpd_buf4_init(&F, CC_FINTS, 0, 28, 26, 28, 26, 0, "F <Ab|Ci>");
    dpd_buf4_init(&L2, CC_LAMPS, 0, 22, 28, 22, 28, 0, "LIjAb");
    dpd_contract444(&L2, &F, &Z, 0, 1, -1, 0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&F);
    /** New L(Ij,Ab) <-- Z(Ij,Am) t(m,b) **/
    dpd_buf4_init(&newLIjAb, CC_LAMPS, 0, 22, 28, 22, 28, 0, "New LIjAb");
    dpd_contract424(&Z, &tia, &newLIjAb, 3, 0, 0, 1, 1);
    dpd_buf4_close(&newLIjAb);
    dpd_buf4_close(&Z);

    /** Z(Ij,Mb) = - L(Ij,Ef) F(Mb,Ef) **/
    dpd_buf4_init(&Z, CC_TMP1, 0, 22, 24, 22, 24, 0, "Z(Ij,Mb)");
    dpd_buf4_init(&F, CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    dpd_buf4_init(&L2, CC_LAMPS, 0, 22, 28, 22, 28, 0, "LIjAb");
    dpd_contract444(&L2, &F, &Z, 0, 0, -1, 0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&F);
    /** New L(Ij,Ab) <-- t(M,A) Z(Ij,Mb) **/
    dpd_buf4_init(&newLIjAb, CC_LAMPS, 0, 22, 28, 22, 28, 0, "New LIjAb");
    dpd_contract244(&tIA, &Z, &newLIjAb, 0, 2, 1, 1, 1);
    dpd_buf4_close(&newLIjAb);
    dpd_buf4_close(&Z);

    /** Z(Ij,Mn) = L(Ij,Ef) <Mn|Ef> **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 22, 22, 22, 22, 0, "Z(Ij,Mn)");
    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_buf4_init(&L2, CC_LAMPS, 0, 22, 28, 22, 28, 0, "LIjAb");
    dpd_contract444(&L2, &D, &Z, 0, 0, 1, 0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&D);
    /** New LIjAb <-- Z(Ij,Mn) Tau(Mn,Ab) **/
    dpd_buf4_init(&newLIjAb, CC_LAMPS, 0, 22, 28, 22, 28, 0, "New LIjAb");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
    dpd_contract444(&Z, &T2, &newLIjAb, 0, 1, 1, 1);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&newLIjAb);
    dpd_buf4_close(&Z);


    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);

  }

}
