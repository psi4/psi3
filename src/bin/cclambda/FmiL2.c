#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/** The RHF/ROHF contractions can be improved here **/

void FmiL2(void)
{
  dpdbuf4 Lijab, LIJAB, LIjAb;
  dpdbuf4 newLijab, newLIJAB, newLIjAb;
  dpdfile2 LFmit2, LFMIt2;
  dpdbuf4 X, X1, X2;

  /* RHS -= P(ij)*Limab*Fjm */
  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/

    dpd_file2_init(&LFMIt2, CC_OEI, 0, 0, 0, "FMIt");
    dpd_file2_init(&LFmit2, CC_OEI, 0, 0, 0, "Fmit");

    dpd_buf4_init(&LIJAB, CC_LAMPS, 0, 0, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&X1, CC_TMP1, 0, 0, 7, 0, 7, 0, "X(0,7) 1");
    dpd_contract424(&LIJAB, &LFMIt2, &X1, 1, 1, 1, -1.0, 0.0);
    dpd_buf4_init(&X2, CC_TMP1, 0, 0, 7, 0, 7, 0, "X(0,7) 2");
    dpd_contract244(&LFMIt2, &LIJAB, &X2, 1, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&LIJAB);
    dpd_buf4_axpy(&X1, &X2, 1.0);
    dpd_buf4_close(&X1);
    dpd_buf4_init(&newLIJAB, CC_LAMPS, 0, 0, 7, 2, 7, 0, "New LIJAB");
    dpd_buf4_axpy(&X2, &newLIJAB, 1.0);
    dpd_buf4_close(&X2);
    dpd_buf4_close(&newLIJAB);

    dpd_buf4_init(&Lijab, CC_LAMPS, 0, 0, 7, 2, 7, 0, "Lijab");
    dpd_buf4_init(&X1, CC_TMP1, 0, 0, 7, 0, 7, 0, "X(0,7) 1");
    dpd_contract424(&Lijab, &LFmit2, &X1, 1, 1, 1, -1.0, 0.0);
    dpd_buf4_init(&X2, CC_TMP1, 0, 0, 7, 0, 7, 0, "X(0,7) 2");
    dpd_contract244(&LFmit2, &Lijab, &X2, 1, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&Lijab);
    dpd_buf4_axpy(&X1, &X2, 1.0);
    dpd_buf4_close(&X1);
    dpd_buf4_init(&newLijab, CC_LAMPS, 0, 0, 7, 2, 7, 0, "New Lijab");
    dpd_buf4_axpy(&X2, &newLijab, 1.0);
    dpd_buf4_close(&X2);
    dpd_buf4_close(&newLijab);

    dpd_buf4_init(&LIjAb, CC_LAMPS, 0, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_init(&newLIjAb, CC_LAMPS, 0, 0, 5, 0, 5, 0, "New LIjAb");
    dpd_contract424(&LIjAb, &LFmit2, &newLIjAb, 1, 1, 1, -1.0, 1.0);
    dpd_contract244(&LFMIt2, &LIjAb, &newLIjAb, 1, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&LIjAb);
    dpd_buf4_close(&newLIjAb);

    dpd_file2_close(&LFmit2);
    dpd_file2_close(&LFMIt2);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_file2_init(&LFMIt2, CC_OEI, 0, 0, 0, "FMIt");
    dpd_file2_init(&LFmit2, CC_OEI, 0, 2, 2, "Fmit");

    /** X(IJ,AB) = F(I,M) L(MJ,AB) **/
    dpd_buf4_init(&X, CC_TMP1, 0, 0, 7, 0, 7, 0, "X(IJ,AB)");
    dpd_buf4_init(&LIJAB, CC_LAMPS, 0, 0, 7, 2, 7, 0, "LIJAB");
    dpd_contract244(&LFMIt2, &LIJAB, &X, 1, 0, 0, -1, 0);
    dpd_buf4_close(&LIJAB);
    /** X(IJ,AB) --> X'(JI,AB) **/
    dpd_buf4_sort(&X, CC_TMP1, qprs, 0, 7, "X'(JI,AB)");
    dpd_buf4_close(&X);

    /** X(IJ,AB) = X(IJ,AB) - X'(JI,AB) **/
    dpd_buf4_init(&X1, CC_TMP1, 0, 0, 7, 0, 7, 0, "X(IJ,AB)");
    dpd_buf4_init(&X2, CC_TMP1, 0, 0, 7, 0, 7, 0, "X'(JI,AB)");
    dpd_buf4_axpy(&X2, &X1, -1.0);
    dpd_buf4_close(&X2);
    /** L(IJ,AB) <--- X(IJ,AB) **/
    dpd_buf4_init(&newLIJAB, CC_LAMPS, 0, 0, 7, 2, 7, 0, "New LIJAB");
    dpd_buf4_axpy(&X1, &newLIJAB, 1.0);
    dpd_buf4_close(&X1);
    dpd_buf4_close(&newLIJAB);


    /** X(ij,ab) = F(i,m) L(mj,ab) **/
    dpd_buf4_init(&X, CC_TMP1, 0, 10, 17, 10, 17, 0, "X(ij,ab)");
    dpd_buf4_init(&LIJAB, CC_LAMPS, 0, 10, 17, 12, 17, 0, "Lijab");
    dpd_contract244(&LFmit2, &LIJAB, &X, 1, 0, 0, -1, 0);
    dpd_buf4_close(&LIJAB);
    /** X(ij,ab) --> X'(ji,ab) **/
    dpd_buf4_sort(&X, CC_TMP1, qprs, 10, 17, "X'(ji,ab)");
    dpd_buf4_close(&X);

    /** X(ij,ab) = X(ij,ab) - X'(ji,ab) **/
    dpd_buf4_init(&X1, CC_TMP1, 0, 10, 17, 10, 17, 0, "X(ij,ab)");
    dpd_buf4_init(&X2, CC_TMP1, 0, 10, 17, 10, 17, 0, "X'(ji,ab)");
    dpd_buf4_axpy(&X2, &X1, -1.0);
    dpd_buf4_close(&X2);
    /** L(ij,ab) <--- X(ij,ab) **/
    dpd_buf4_init(&newLIJAB, CC_LAMPS, 0, 10, 17, 12, 17, 0, "New Lijab");
    dpd_buf4_axpy(&X1, &newLIJAB, 1.0);
    dpd_buf4_close(&X1);
    dpd_buf4_close(&newLIJAB);

    /** L(Ij,Ab) <-- L(Im,Ab) F(j,m) - F(I,M) L(Mj,Ab) **/
    dpd_buf4_init(&LIjAb, CC_LAMPS, 0, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_init(&newLIjAb, CC_LAMPS, 0, 22, 28, 22, 28, 0, "New LIjAb");
    dpd_contract424(&LIjAb, &LFmit2, &newLIjAb, 1, 1, 1, -1, 1);
    dpd_contract244(&LFMIt2, &LIjAb, &newLIjAb, 1, 0, 0, -1, 1);
    dpd_buf4_close(&LIjAb);
    dpd_buf4_close(&newLIjAb);

    dpd_file2_close(&LFmit2);
    dpd_file2_close(&LFMIt2);
  }
}
