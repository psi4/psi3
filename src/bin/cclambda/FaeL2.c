#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/** The ROHF version of this contraction can be done with fewer contractions.
 **/

void FaeL2(void)
{
  dpdbuf4 Lijab, LIJAB, LIjAb;
  dpdbuf4 newLIJAB, newLijab, newLIjAb;
  dpdfile2 LFaet2, LFAEt2;
  dpdbuf4 X, X1, X2;

  /* RHS += P(ab)*Lijae*Feb */

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/

    dpd_file2_init(&LFAEt2, CC_OEI, 0, 1, 1, "FAEt");
    dpd_file2_init(&LFaet2, CC_OEI, 0, 1, 1, "Faet");

    dpd_buf4_init(&LIJAB, CC_LAMPS, 0, 2, 5, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&X1, CC_TMP1, 0, 2, 5, 2, 5, 0, "X(2,5) 1");
    dpd_contract424(&LIJAB, &LFAEt2, &X1, 3, 0, 0, 1.0, 0.0);
    dpd_buf4_init(&X2, CC_TMP1, 0, 2, 5, 2, 5, 0, "X(2,5) 2");
    dpd_contract244(&LFAEt2, &LIJAB, &X2, 0, 2, 1, 1.0, 0.0);
    dpd_buf4_close(&LIJAB);
    dpd_buf4_axpy(&X1, &X2, 1.0);
    dpd_buf4_close(&X1);
    dpd_buf4_init(&newLIJAB, CC_LAMPS, 0, 2, 5, 2, 7, 0, "New LIJAB");
    dpd_buf4_axpy(&X2, &newLIJAB, 1.0);
    dpd_buf4_close(&X2);
    dpd_buf4_close(&newLIJAB);

    dpd_buf4_init(&Lijab, CC_LAMPS, 0, 2, 5, 2, 7, 0, "Lijab");
    dpd_buf4_init(&X1, CC_TMP1, 0, 2, 5, 2, 5, 0, "X(2,5) 1");
    dpd_contract424(&Lijab, &LFaet2, &X1, 3, 0, 0, 1.0, 0.0);
    dpd_buf4_init(&X2, CC_TMP1, 0, 2, 5, 2, 5, 0, "X(2,5) 2");
    dpd_contract244(&LFaet2, &Lijab, &X2, 0, 2, 1, 1.0, 0.0);
    dpd_buf4_close(&Lijab);
    dpd_buf4_axpy(&X1, &X2, 1.0);
    dpd_buf4_close(&X1);
    dpd_buf4_init(&newLijab, CC_LAMPS, 0, 2, 5, 2, 7, 0, "New Lijab");
    dpd_buf4_axpy(&X2, &newLijab, 1.0);
    dpd_buf4_close(&X2);
    dpd_buf4_close(&newLijab);

    dpd_buf4_init(&LIjAb, CC_LAMPS, 0, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_init(&newLIjAb, CC_LAMPS, 0, 0, 5, 0, 5, 0, "New LIjAb");
    dpd_contract424(&LIjAb, &LFaet2, &newLIjAb, 3, 0, 0, 1.0, 1.0);
    dpd_contract244(&LFAEt2, &LIjAb, &newLIjAb, 0, 2, 1, 1.0, 1.0);
    dpd_buf4_close(&LIjAb);
    dpd_buf4_close(&newLIjAb);

    dpd_file2_close(&LFaet2);
    dpd_file2_close(&LFAEt2);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_file2_init(&LFAEt2, CC_OEI, 0, 1, 1, "FAEt");
    dpd_file2_init(&LFaet2, CC_OEI, 0, 3, 3, "Faet");

    /** X(IJ,AB) = L_IJ^AE F_EB **/
    dpd_buf4_init(&X, CC_TMP1, 0, 2, 5, 2, 5, 0, "X(IJ,AB)");
    dpd_buf4_init(&LIJAB, CC_LAMPS, 0, 2, 5, 2, 7, 0, "LIJAB");
    dpd_contract424(&LIJAB, &LFAEt2, &X, 3, 0, 0, 1, 0);
    dpd_buf4_close(&LIJAB);
    /** X(IJ,AB) --> X'(IJ,BA) **/
    dpd_buf4_sort(&X, CC_TMP1, pqsr, 2, 5, "X'(IJ,BA)");
    dpd_buf4_close(&X);
    /** X(IJ,AB) = X(IJ,AB) - X'(IJ,BA) **/
    dpd_buf4_init(&X1, CC_TMP1, 0, 2, 5, 2, 5, 0, "X(IJ,AB)");
    dpd_buf4_init(&X2, CC_TMP1, 0, 2, 5, 2, 5, 0, "X'(IJ,BA)");
    dpd_buf4_axpy(&X2, &X1, -1);
    dpd_buf4_close(&X2);
    dpd_buf4_close(&X1);
    /** L(IJ,AB) <-- X(IJ,AB) **/
    dpd_buf4_init(&X, CC_TMP1, 0, 2, 5, 2, 5, 0, "X(IJ,AB)");
    dpd_buf4_init(&newLIJAB, CC_LAMPS, 0, 2, 5, 2, 7, 0, "New LIJAB");
    dpd_buf4_axpy(&X, &newLIJAB, 1.0);
    dpd_buf4_close(&X);
    dpd_buf4_close(&newLIJAB);

    /** X(ij,ab) = L_ij^ae F_eb **/
    dpd_buf4_init(&X, CC_TMP1, 0, 12, 15, 12, 15, 0, "X(ij,ab)");
    dpd_buf4_init(&LIJAB, CC_LAMPS, 0, 12, 15, 12, 17, 0, "Lijab");
    dpd_contract424(&LIJAB, &LFaet2, &X, 3, 0, 0, 1, 0);
    dpd_buf4_close(&LIJAB);
    /** X(ij,ab) --> X'(ij,ba) **/
    dpd_buf4_sort(&X, CC_TMP1, pqsr, 12, 15, "X'(ij,ba)");
    dpd_buf4_close(&X);
    /** X(ij,ab) = X(ij,ab) - X'(ij,ba) **/
    dpd_buf4_init(&X1, CC_TMP1, 0, 12, 15, 12, 15, 0, "X(ij,ab)");
    dpd_buf4_init(&X2, CC_TMP1, 0, 12, 15, 12, 15, 0, "X'(ij,ba)");
    dpd_buf4_axpy(&X2, &X1, -1);
    dpd_buf4_close(&X2);
    dpd_buf4_close(&X1);
    /** L(ij,ab) <-- X(ij,ab) **/
    dpd_buf4_init(&X, CC_TMP1, 0, 12, 15, 12, 15, 0, "X(ij,ab)");
    dpd_buf4_init(&newLIJAB, CC_LAMPS, 0, 12, 15, 12, 17, 0, "New Lijab");
    dpd_buf4_axpy(&X, &newLIJAB, 1.0);
    dpd_buf4_close(&X);
    dpd_buf4_close(&newLIJAB);

    /** L(Ij,Ab) <-- L(Ij,Ae) F(e,b) - F(E,A) L(Ij,Eb) **/
    dpd_buf4_init(&LIjAb, CC_LAMPS, 0, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_init(&newLIjAb, CC_LAMPS, 0, 22, 28, 22, 28, 0, "New LIjAb");
    dpd_contract424(&LIjAb, &LFaet2, &newLIjAb, 3, 0, 0, 1, 1);
    dpd_contract244(&LFAEt2, &LIjAb, &newLIjAb, 0, 2, 1, 1, 1);
    dpd_buf4_close(&LIjAb);
    dpd_buf4_close(&newLIjAb);

    dpd_file2_close(&LFaet2);
    dpd_file2_close(&LFAEt2);

  }

}
