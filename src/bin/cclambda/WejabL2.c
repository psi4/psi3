#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/* WejabL2(): Computes the contribution of the Wamef HBAR matrix
** elements to the Lambda double de-excitation amplitude equations.
** These contributions are given in spin orbitals as:
**
** L_ij^ab = P(ij) L_i^e Wejab
**
** where Wejab = Wamef and is defined as:
**
** Wamef = <am||ef> - t_n^a <nm||ef>
**
** [cf. Gauss and Stanton, JCP 103, 3561-3577 (1995).]
**
** By spin case of L_ij^ab, these contributions are computed as:
**
** L(IJ,AB) = L(I,E) W(EJ,AB) - L(J,E) W(EI,AB) (only one unique contraction)
** L(ij,ab) = L(i,e) W(ej,ab) - L(j,e) W(ei,ab) (only one unique contraction)
** L(Ij,Ab) = L(I,E) W(Ej,Ab) + L(j,e) W(eI,bA)
**
** TDC, July 2002
**
** NB: The ROHF case needs to be re-written to use (AM,EF) ordering
** for the Wamef matrix elements, as I've done for the UHF case.
*/

void WejabL2(void)
{
  dpdbuf4 W, Wamef, WAmEf, WaMeF, WAMEF;
  dpdbuf4 L2, newLijab, newLIJAB, newLIjAb;
  dpdbuf4 Ltmp;
  dpdfile2 LIA, Lia;
  dpdbuf4 X1, X2, Z, Z1, Z2;
  
  /* RHS += P(ij) Lie * Wejab */
  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/

    dpd_file2_init(&LIA, CC_OEI, L_irr, 0, 1, "LIA");
    dpd_file2_init(&Lia, CC_OEI, L_irr, 0, 1, "Lia");

    dpd_buf4_init(&WAMEF, CC_HBAR, 0, 10, 7, 10, 7, 0, "WAMEF");
    dpd_buf4_init(&X1, CC_TMP1, L_irr, 0, 7, 0, 7, 0, "X(0,7) 1");
    dpd_contract424(&WAMEF, &LIA, &X1, 1, 1, 1, -1.0, 0.0);
    dpd_buf4_close(&WAMEF);
    dpd_buf4_sort(&X1, CC_TMP1, qprs, 0, 7, "X(0,7) 2");
    dpd_buf4_init(&X2, CC_TMP1, L_irr, 0, 7, 0, 7, 0, "X(0,7) 2");
    dpd_buf4_axpy(&X2, &X1, -1.0);
    dpd_buf4_close(&X2);
    dpd_buf4_init(&newLIJAB, CC_LAMPS, L_irr, 0, 7, 2, 7, 0, "New LIJAB");
    dpd_buf4_axpy(&X1, &newLIJAB, 1.0);
    dpd_buf4_close(&X1);
    dpd_buf4_close(&newLIJAB);

    /*
      dpd_buf4_init(&WAMEF, CC_HBAR, 0, 10, 7, 10, 7, 0, "WAMEF");
      dpd_buf4_init(&newLIJAB, CC_LAMPS, L_irr, 0, 7, 2, 7, 0, "New LIJAB");
      dpd_contract424(&WAMEF, &LIA, &newLIJAB, 1, 1, 1, -1.0, 1.0);

      dpd_buf4_init(&Ltmp, CC_TMP0, L_irr, 0, 7, 0, 7, 0, "LIJAB (JI,A>B)");
      dpd_contract424(&WAMEF, &LIA, &Ltmp, 1, 1, 1, 1.0, 0.0);
      dpd_buf4_sort(&Ltmp, CC_TMP1, qprs, 0, 7, "LIJAB (IJ,A>B)");
      dpd_buf4_close(&Ltmp);
      dpd_buf4_init(&Ltmp, CC_TMP1, L_irr, 0, 7, 0, 7, 0, "LIJAB (IJ,A>B)");
      dpd_buf4_axpy(&Ltmp, &newLIJAB, 1.0);
      dpd_buf4_close(&Ltmp);
      dpd_buf4_close(&newLIJAB);
      dpd_buf4_close(&WAMEF);
    */

    dpd_buf4_init(&Wamef, CC_HBAR, 0, 10, 7, 10, 7, 0, "Wamef");
    dpd_buf4_init(&X1, CC_TMP1, L_irr, 0, 7, 0, 7, 0, "X(0,7) 1");
    dpd_contract424(&Wamef, &Lia, &X1, 1, 1, 1, -1.0, 0.0);
    dpd_buf4_close(&Wamef);
    dpd_buf4_sort(&X1, CC_TMP1, qprs, 0, 7, "X(0,7) 2");
    dpd_buf4_init(&X2, CC_TMP1, L_irr, 0, 7, 0, 7, 0, "X(0,7) 2");
    dpd_buf4_axpy(&X2, &X1, -1.0);
    dpd_buf4_close(&X2);
    dpd_buf4_init(&newLijab, CC_LAMPS, L_irr, 0, 7, 2, 7, 0, "New Lijab");
    dpd_buf4_axpy(&X1, &newLijab, 1.0);
    dpd_buf4_close(&X1);
    dpd_buf4_close(&newLijab);

    /*
      dpd_buf4_init(&Wamef, CC_HBAR, 0, 10, 7, 10, 7, 0, "Wamef");
      dpd_buf4_init(&newLijab, CC_LAMPS, L_irr, 0, 7, 2, 7, 0, "New Lijab");
      dpd_contract424(&Wamef, &Lia, &newLijab, 1, 1, 1, -1.0, 1.0);

      dpd_buf4_init(&Ltmp, CC_TMP0, L_irr, 0, 7, 0, 7, 0, "Lijab (ji,a>b)");
      dpd_contract424(&Wamef, &Lia, &Ltmp, 1, 1, 1, 1.0, 0.0);
      dpd_buf4_sort(&Ltmp, CC_TMP1, qprs, 0, 7, "Lijab (ij,a>b)");
      dpd_buf4_close(&Ltmp);
      dpd_buf4_init(&Ltmp, CC_TMP1, L_irr, 0, 7, 0, 7, 0, "Lijab (ij,a>b)");
      dpd_buf4_axpy(&Ltmp, &newLijab, 1.0);
      dpd_buf4_close(&Ltmp);
      dpd_buf4_close(&newLijab);
      dpd_buf4_close(&Wamef);
    */

    dpd_buf4_init(&newLIjAb, CC_LAMPS, L_irr, 0, 5, 0, 5, 0, "New LIjAb");

    dpd_buf4_init(&WaMeF, CC_HBAR, 0, 10, 5, 10, 5, 0, "WaMeF");
    dpd_buf4_sort(&WaMeF, CC_TMP0, pqsr, 10, 5, "WaMeF (Ma,Fe)");
    dpd_buf4_close(&WaMeF);

    dpd_buf4_init(&WaMeF, CC_TMP0, 0, 10, 5, 10, 5, 0, "WaMeF (Ma,Fe)");
    dpd_contract424(&WaMeF, &Lia, &newLIjAb, 1, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&WaMeF);

    dpd_buf4_init(&WAmEf, CC_HBAR, 0, 10, 5, 10, 5, 0, "WAmEf");
    dpd_buf4_init(&Ltmp, CC_TMP0, L_irr, 0, 5, 0, 5, 0, "Ltmp (ji,ab)");
    dpd_contract424(&WAmEf, &LIA, &Ltmp, 1, 1, 1, 1.0, 0.0);
    dpd_buf4_sort(&Ltmp, CC_TMP1, qprs, 0, 5, "Lijab (ij,ab)");
    dpd_buf4_close(&Ltmp);
    dpd_buf4_close(&WAmEf);

    dpd_buf4_init(&Ltmp, CC_TMP1, L_irr, 0, 5, 0, 5, 0, "Lijab (ij,ab)");
    dpd_buf4_axpy(&Ltmp, &newLIjAb, 1.0);
    dpd_buf4_close(&Ltmp);

    dpd_buf4_close(&newLIjAb);
    dpd_file2_close(&Lia);
    dpd_file2_close(&LIA);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_file2_init(&LIA, CC_OEI, L_irr, 0, 1, "LIA");
    dpd_file2_init(&Lia, CC_OEI, L_irr, 2, 3, "Lia");

    /** Z(IJ,AB) = L(I,E) W(EJ,AB) **/
    dpd_buf4_init(&Z, CC_TMP1, L_irr, 0, 7, 0, 7, 0, "Z(IJ,AB)");
    dpd_buf4_init(&W, CC_HBAR, 0, 21, 7, 21, 7, 0, "WAMEF");
    dpd_contract244(&LIA, &W, &Z, 1, 0, 0, 1, 0);
    dpd_buf4_close(&W);
    /** Z(IJ,AB) --> Z(JI,AB) **/
    dpd_buf4_sort(&Z, CC_TMP1, qprs, 0, 7, "Z(JI,AB)");
    dpd_buf4_close(&Z);
    /** Z(IJ,AB) = Z(IJ,AB) - Z(JI,AB) **/
    dpd_buf4_init(&Z1, CC_TMP1, L_irr, 0, 7, 0, 7, 0, "Z(IJ,AB)");
    dpd_buf4_init(&Z2, CC_TMP1, L_irr, 0, 7, 0, 7, 0, "Z(JI,AB)");
    dpd_buf4_axpy(&Z2, &Z1, -1);
    dpd_buf4_close(&Z2);
    /** Z(IJ,AB) --> New L(IJ,AB) **/
    dpd_buf4_init(&L2, CC_LAMPS, L_irr, 0, 7, 2, 7, 0, "New LIJAB");
    dpd_buf4_axpy(&Z1, &L2, 1);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&Z1);

    /** Z(ij,ab) = L(i,e) W(ej,ab) **/
    dpd_buf4_init(&Z, CC_TMP1, L_irr, 10, 17, 10, 17, 0, "Z(ij,ab)");
    dpd_buf4_init(&W, CC_HBAR, 0, 31, 17, 31, 17, 0, "Wamef");
    dpd_contract244(&Lia, &W, &Z, 1, 0, 0, 1, 0);
    dpd_buf4_close(&W);
    /** Z(ij,ab) --> Z(ji,ab) **/
    dpd_buf4_sort(&Z, CC_TMP1, qprs, 10, 17, "Z(ji,ab)");
    dpd_buf4_close(&Z);
    /** Z(ij,ab) = Z(ij,ab) - Z(ji,ab) **/
    dpd_buf4_init(&Z1, CC_TMP1, L_irr, 10, 17, 10, 17, 0, "Z(ij,ab)");
    dpd_buf4_init(&Z2, CC_TMP1, L_irr, 10, 17, 10, 17, 0, "Z(ji,ab)");
    dpd_buf4_axpy(&Z2, &Z1, -1);
    dpd_buf4_close(&Z2);
    /** Z(ij,ab) --> New L(ij,ab) **/
    dpd_buf4_init(&L2, CC_LAMPS, L_irr, 10, 17, 12, 17, 0, "New Lijab");
    dpd_buf4_axpy(&Z1, &L2, 1);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&Z1);


    /** New L(Ij,Ab) <-- L(I,E) W(Ej,Ab) **/
    dpd_buf4_init(&L2, CC_LAMPS, L_irr, 22, 28, 22, 28, 0, "New LIjAb");
    dpd_buf4_init(&W, CC_HBAR, 0, 26, 28, 26, 28, 0, "WAmEf");
    dpd_contract244(&LIA, &W, &L2, 1, 0, 0, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&L2);

    /** Z(jI,bA) = -L(j,e) W(eI,bA) **/
    dpd_buf4_init(&Z, CC_TMP1, L_irr, 23, 29, 23, 29, 0, "Z(jI,bA)");
    dpd_buf4_init(&W, CC_HBAR, 0, 25, 29, 25, 29, 0, "WaMeF");
    dpd_contract244(&Lia, &W, &Z, 1, 0, 0, 1, 0);
    dpd_buf4_close(&W);
    /** Z(jI,bA) --> New L(Ij,Ab) **/
    dpd_buf4_sort_axpy(&Z, CC_LAMPS, qpsr, 22, 28, "New LIjAb", 1);
    dpd_buf4_close(&Z);

    dpd_file2_close(&Lia);
    dpd_file2_close(&LIA);
  }
}
