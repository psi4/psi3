#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/* WefabL2(): Computes the contribution of the Wefab HBAR matrix
** elements to the Lambda double de-excitation amplitude equations.
** These contributions are given in spin-orbitals as:
** 
** L_ij^ab = 1/2 L_ij^ef Wefab
**
** where W_efab is defined in spin orbitals as:
**
** Wefab = <ef||ab> - P(ef) t_m^f <em||ab> + 1/2 tau_mn^ef <mn||ab>
**
** and tau_mn^ef = t_mn^ef + t_m^e t_n^f - t_m^f t_n^e.
**
** [cf. Gauss and Stanton, JCP 103, 3561-3577 (1995).]
**
** NB: Wefab is not symmetric, so one must be careful when defining
** intermediate quantities for efficient contractions.  I use the
** following contraction steps for each spin case:
**
** Wefab term II spin cases: 
**
**   L_IJ^AB <-- 1/2 ( -t_M^F <EM||AB> L_IJ^EF + t_M^E <FM||AB> L_IJ^EF )
**
**     Z(IJ,EM) = -t_M^F L_IJ^EF
**
**   L_IJ^AB <-- Z(IJ,EM) <EM||AB>
*******
**   L_ij^ab <-- 1/2 ( -t_m^f <em||ab> L_ij^ef + t_m^e <fm||ab> L_ij^ef )
**
**     Z(ij,em) = -t_m^f L_ij^ef
**
**   L_ij^ab <-- Z(ij,em) <em||ab>
*******
**   L_Ij^Ab <-- -t_m^f <Em|Ab> L_Ij^Ef - t_M^E <Mf|Ab>  L_Ij^Ef
**
**     Z(Ij,Em) = -t_m^f L_Ij^Ef
**     Z(Ij,Mf) = -t_M^E L_Ij^Ef
**
**   L_Ij^Ab <-- Z(Ij,Em) <Em|Ab> + Z(Ij,Mf) <Mf|Ab>
**
** Wefab term III:
**
**   L_IJ^AB <-- 1/4 tau_MN^EF <MN||AB> L_IJ^EF
**
**     Z(IJ,MN) = 1/2 tau_MN^EF L_IJ^EF
**
**   L_IJ^AB <-- 1/2 Z(IJ,MN) <MN||AB>
*******
**   L_ij^ab <-- 1/4 tau_mn^ef <mn||ab> L_ij^ef
**
**     Z(ij,mn) = 1/2 tau_mn^ef L_ij^ef
**
**   L_ij^ab <-- 1/2 Z(ij,mn) <mn||ab>
*******
**   L_Ij^Ab <-- tau_Mn^Ef <Mn|Ab> L_Ij^Ef
**
**     Z(Ij,Mn) = tau_Mn^Ef L_Ij^Ef
**
**   L_Ij^Ab <-- Z(Ij,Mn) <Mn|Ab>
*******
**
** TDC, July 2002
*/

void WefabL2(int L_irr)
{
  dpdbuf4 Lijab, LIJAB, LIjAb;
  dpdbuf4 newLijab, newLIJAB, newLIjAb;
  dpdbuf4 Tau, T2, Z, Z1, Z2, L, L2, B, D, F, Ltmp;
  dpdfile2 tIA, tia;

  /* RHS += Wefab*Lijef  */
  if(params.ref == 0) { /** RHF **/

    dpd_buf4_init(&LIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_init(&Z, CC_TMP0, L_irr, 5, 0, 5, 0, 0, "ZAbIj");
    dpd_buf4_init(&B, CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
    dpd_contract444(&B, &LIjAb, &Z, 0, 0, 1, 0);
    dpd_buf4_close(&B);
    dpd_buf4_sort_axpy(&Z, CC_LAMBDA, rspq, 0, 5, "New LIjAb", 1);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&LIjAb);

    dpd_buf4_init(&newLIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_buf4_init(&LIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_init(&Z, CC_TMP0, L_irr, 10, 0, 10, 0, 0, "Z(Mf,Ij)");
    dpd_contract244(&tIA, &LIjAb, &Z, 1, 2, 0, 1, 0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&LIjAb);
    dpd_file2_close(&tIA);

    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_buf4_init(&Z, CC_TMP0, L_irr, 10, 0, 10, 0, 0, "Z(Mf,Ij)");
    dpd_buf4_init(&Z1, CC_TMP0, L_irr, 5, 0, 5, 0, 0, "Z(Ab,Ij)");
    dpd_contract444(&F, &Z, &Z1, 1, 1, -1, 0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&F);
    dpd_buf4_close(&newLIjAb);
    dpd_buf4_sort_axpy(&Z1, CC_LAMBDA, srqp, 0, 5, "New LIjAb", 1);
    dpd_buf4_sort_axpy(&Z1, CC_LAMBDA, rspq, 0, 5, "New LIjAb", 1);
    dpd_buf4_init(&newLIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    dpd_buf4_close(&Z1);

    dpd_buf4_init(&Z, CC_TMP0, L_irr, 0, 0, 0, 0, 0, "Z(Ij,Mn)");
    dpd_buf4_init(&Tau, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_contract444(&L2, &Tau, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&Tau);
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract444(&Z, &D, &newLIjAb, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&Z);

    dpd_buf4_close(&newLIjAb);

  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    dpd_buf4_init(&LIJAB, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&Z, CC_TMP2, L_irr, 7, 2, 7, 2, 0, "ZABIJ");
    dpd_buf4_init(&B, CC_BINTS, 0, 7, 7, 5, 5, 1, "B <ab|cd>");
    dpd_contract444(&B, &LIJAB, &Z, 0, 0, 1, 0);
    dpd_buf4_close(&B);
    dpd_buf4_sort_axpy(&Z, CC_LAMBDA, rspq, 2, 7, "New LIJAB", 1);
    dpd_buf4_close(&LIJAB);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&newLIJAB, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");

    dpd_buf4_init(&LIJAB, CC_LAMBDA, L_irr, 2, 5, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&Ltmp, CC_TMP0, L_irr, 2, 10, 2, 10, 0, "Ltmp (I>J,MF)");
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
    dpd_contract244(&tIA, &LIJAB, &Ltmp, 1, 2, 1, 1.0, 0.0);
    dpd_contract444(&Ltmp, &F, &newLIJAB, 0, 1, -1.0, 1.0);
    dpd_buf4_close(&F);
    dpd_buf4_close(&Ltmp);
    dpd_buf4_close(&LIJAB);

    dpd_buf4_init(&Z, CC_TMP0, L_irr, 2, 2, 2, 2, 0, "Z(IJ,MN)");
    dpd_buf4_init(&Tau, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    dpd_contract444(&L2, &Tau, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&Tau);
    dpd_buf4_init(&D, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
    dpd_contract444(&Z, &D, &newLIJAB, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&newLIJAB);

    dpd_buf4_init(&Lijab, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
    dpd_buf4_init(&Z, CC_TMP2, L_irr, 7, 2, 7, 2, 0, "Zabij");
    dpd_buf4_init(&B, CC_BINTS, 0, 7, 7, 5, 5, 1, "B <ab|cd>");
    dpd_contract444(&B, &Lijab, &Z, 0, 0, 1, 0);
    dpd_buf4_close(&B);
    dpd_buf4_close(&Lijab);
    dpd_buf4_sort_axpy(&Z, CC_LAMBDA, rspq, 2, 7, "New Lijab", 1);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&newLijab, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab");

    dpd_buf4_init(&Lijab, CC_LAMBDA, L_irr, 2, 5, 2, 7, 0, "Lijab");
    dpd_buf4_init(&Ltmp, CC_TMP0, L_irr, 2, 10, 2, 10, 0, "Ltmp (i>j,mf)");
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
    dpd_contract244(&tia, &Lijab, &Ltmp, 1, 2, 1, 1.0, 0.0);
    dpd_contract444(&Ltmp, &F, &newLijab, 0, 1, -1.0, 1.0);
    dpd_buf4_close(&F);
    dpd_buf4_close(&Ltmp);
    dpd_buf4_close(&Lijab);

    dpd_buf4_init(&Z, CC_TMP0, L_irr, 2, 2, 2, 2, 0, "Z(ij,mn)");
    dpd_buf4_init(&Tau, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
    dpd_contract444(&L2, &Tau, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&Tau);
    dpd_buf4_init(&D, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
    dpd_contract444(&Z, &D, &newLijab, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&newLijab);


    dpd_buf4_init(&LIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_init(&Z, CC_TMP2, L_irr, 5, 0, 5, 0, 0, "ZAbIj");
    dpd_buf4_init(&B, CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
    dpd_contract444(&B, &LIjAb, &Z, 0, 0, 1, 0);
    dpd_buf4_close(&B);
    dpd_buf4_sort_axpy(&Z, CC_LAMBDA, rspq, 0, 5, "New LIjAb", 1);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&newLIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");

    dpd_buf4_init(&Ltmp, CC_TMP1, L_irr, 0, 11, 0, 11, 0, "Lt (Ij,Em)");
    dpd_contract424(&LIjAb, &tia, &Ltmp, 3, 1, 0, 1.0, 0.0);
    dpd_buf4_sort(&Ltmp, CC_TMP2, pqsr, 0, 10, "Lt (Ij,mE)");
    dpd_buf4_close(&Ltmp);
    dpd_buf4_init(&Ltmp, CC_TMP3, L_irr, 0, 10, 0, 10, 0, "Lt (Ij,Mf)");
    dpd_contract244(&tIA, &LIjAb, &Ltmp, 1, 2, 1, 1.0, 0.0);
    dpd_buf4_close(&Ltmp);

    dpd_buf4_close(&LIjAb);

    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_buf4_init(&Ltmp, CC_TMP3, L_irr, 0, 10, 0, 10, 0, "Lt (Ij,Mf)");
    dpd_contract444(&Ltmp, &F, &newLIjAb, 0, 1, -1.0, 1.0);
    dpd_buf4_close(&Ltmp);
    dpd_buf4_sort(&F, CC_TMP0, pqsr, 10, 5, "<me|ab> (mE,Ab)");
    dpd_buf4_close(&F);

    dpd_buf4_init(&F, CC_TMP0, 0, 10, 5, 10, 5, 0, "<me|ab> (mE,Ab)");
    dpd_buf4_init(&Ltmp, CC_TMP2, L_irr, 0, 10, 0, 10, 0, "Lt (Ij,mE)");
    dpd_contract444(&Ltmp, &F, &newLIjAb, 0, 1, -1.0, 1.0);
    dpd_buf4_close(&Ltmp);
    dpd_buf4_close(&F);

    dpd_buf4_init(&Z, CC_TMP0, L_irr, 0, 0, 0, 0, 0, "Z(Ij,Mn)");
    dpd_buf4_init(&Tau, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
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
    dpd_buf4_init(&Z, CC_TMP1, L_irr, 7, 2, 7, 2, 0, "Z(AB,IJ)");
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&B, CC_BINTS, 0, 7, 7, 5, 5, 1, "B <AB|CD>");
    dpd_contract444(&B, &L2, &Z, 0, 0, 1, 0);
    dpd_buf4_close(&B);
    dpd_buf4_close(&L2);
    /** Z(AB,IJ) --> New L(IJ,AB) **/
    dpd_buf4_sort_axpy(&Z, CC_LAMBDA, rspq, 2, 7, "New LIJAB", 1);
    dpd_buf4_close(&Z);

    /** Z(IJ,EM) = -L(IJ,EFf) t(M,F) **/
    dpd_buf4_init(&Z, CC_TMP1, L_irr, 2, 21, 2, 21, 0, "Z(IJ,EM)");
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 2, 5, 2, 7, 0, "LIJAB");
    dpd_contract424(&L2, &tIA, &Z, 3, 1, 0, -1, 0);
    dpd_buf4_close(&L2);
    /** New L(IJ,AB) <-- Z(IJ,EM) <EM||AB> **/
    dpd_buf4_init(&F, CC_FINTS, 0, 21, 7, 21, 5, 1, "F <AI|BC>");
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    dpd_contract444(&Z, &F, &L2, 0, 1, 1, 1);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&F);
    dpd_buf4_close(&Z);

    /** Z(IJ,MN) = 1/2 L(IJ,EF) tau_MN^EF **/
    dpd_buf4_init(&Z, CC_TMP1, L_irr, 2, 2, 2, 2, 0, "Z(IJ,MN)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    dpd_contract444(&L2, &T2, &Z, 0, 0, 1, 0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    /** New L(IJ,AB) <-- 1/2 Z(IJ,MN) <MN||AB> **/
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    dpd_buf4_init(&D, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
    dpd_contract444(&Z, &D, &L2, 0, 1, 1, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&Z);


    /** Z(ab,ij) = <ab||cd> L(ij,cd) **/
    dpd_buf4_init(&Z, CC_TMP1, L_irr, 17, 12, 17, 12, 0, "Z(ab,ij)");
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "Lijab");
    dpd_buf4_init(&B, CC_BINTS, 0, 17, 17, 15, 15, 1, "B <ab|cd>");
    dpd_contract444(&B, &L2, &Z, 0, 0, 1, 0);
    dpd_buf4_close(&B);
    dpd_buf4_close(&L2);
    /** Z(ab,ij) --> New L(ij,ab) **/
    dpd_buf4_sort_axpy(&Z, CC_LAMBDA, rspq, 12, 17, "New Lijab", 1);
    dpd_buf4_close(&Z);

    /** Z(ij,em) = -L(ij,ef) t(m,f) **/
    dpd_buf4_init(&Z, CC_TMP1, L_irr, 12, 31, 12, 31, 0, "Z(ij,em)");
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 12, 15, 12, 17, 0, "Lijab");
    dpd_contract424(&L2, &tia, &Z, 3, 1, 0, -1, 0);
    dpd_buf4_close(&L2);
    /** New L(ij,ab) <-- Z(ij,em) <em||ab> **/
    dpd_buf4_init(&F, CC_FINTS, 0, 31, 17, 31, 15, 1, "F <ai|bc>");
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "New Lijab");
    dpd_contract444(&Z, &F, &L2, 0, 1, 1, 1);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&F);
    dpd_buf4_close(&Z);

    /** Z(ij,mn) = 1/2 L(ij,ef) tau_mn^ef **/
    dpd_buf4_init(&Z, CC_TMP1, L_irr, 12, 12, 12, 12, 0, "Z(ij,mn)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "Lijab");
    dpd_contract444(&L2, &T2, &Z, 0, 0, 1, 0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    /** New L(ij,ab) <-- 1/2 Z(ij,mn) <mn||ab> **/
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "New Lijab");
    dpd_buf4_init(&D, CC_DINTS, 0, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
    dpd_contract444(&Z, &D, &L2, 0, 1, 1, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&Z);


    /** Z(Ab,Ij) = <Ab|Cd> L(Ij,Cd) **/
    dpd_buf4_init(&Z, CC_TMP1, L_irr, 28, 22, 28, 22, 0, "Z(Ab,Ij)");
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_init(&B, CC_BINTS, 0, 28, 28, 28, 28, 0, "B <Ab|Cd>");
    dpd_contract444(&B, &L2, &Z, 0, 0, 1, 0);
    dpd_buf4_close(&B);
    dpd_buf4_close(&L2);
    /** Z(Ab,Ij) --> New L(Ij,Ab) **/
    dpd_buf4_sort_axpy(&Z, CC_LAMBDA, rspq, 22, 28, "New LIjAb", 1);
    dpd_buf4_close(&Z);

    /** Z(Ij,Em) = -L(Ij,Ef) t(m,f) **/
    dpd_buf4_init(&Z, CC_TMP1, L_irr, 22, 26, 22, 26, 0, "Z(Ij,Em)");
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    dpd_contract424(&L2, &tia, &Z, 3, 1, 0, -1, 0);
    dpd_buf4_close(&L2);
    /** New L(Ij,Ab) <-- Z(Ij,Em) <Em|Ab> **/
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "New LIjAb");
    dpd_buf4_init(&F, CC_FINTS, 0, 26, 28, 26, 28, 0, "F <Ai|Bc>");
    dpd_contract444(&Z, &F, &L2, 0, 1, 1, 1);
    dpd_buf4_close(&F);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&Z);

    /** Z(Ij,Mf) = -t(M,E) L(Ij,Ef) **/
    dpd_buf4_init(&Z, CC_TMP1, L_irr, 22, 24, 22, 24, 0, "Z(Ij,Mf)");
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    dpd_contract244(&tIA, &L2, &Z, 1, 2, 1, -1, 0);
    dpd_buf4_close(&L2);
    /** New L(Ij,Ab) <-- Z(Ij,Mf) <Mf|Ab> **/
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "New LIjAb");
    dpd_buf4_init(&F, CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    dpd_contract444(&Z, &F, &L2, 0, 1, 1, 1);
    dpd_buf4_close(&F);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&Z);

    /** Z(Ij,Mn) = L(Ij,Ef) tau(Mn,Ef) **/
    dpd_buf4_init(&Z, CC_TMP1, L_irr, 22, 22, 22, 22, 0, "Z(Ij,Mn)");
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
    dpd_contract444(&L2, &T2, &Z, 0, 0, 1, 0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    /** New L(Ij,Ab) <-- Z(Ij,Mn) <Mn|Ab> **/
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "New LIjAb");
    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_contract444(&Z, &D, &L2, 0, 1, 1, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&Z);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);

  }

}
