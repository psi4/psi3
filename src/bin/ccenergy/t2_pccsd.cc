/*! \file
    \ingroup CCENERGY
    \brief implements extra terms from pCCSD theory of Nooijen
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <assert.h>
#include <libdpd/dpd.h>
#include "Params.h"
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

void status(const char *, FILE *);

void t2_delta_pCCSD_build(void)
{
  // only RHF implemented
  assert( params.ref == 0 );

  if (params.pccsd_alpha != 1.0) {

    //
    // A = v^{kl}_{cd}*  t2^{cd}_{jl} * t2^{ab}_{ik} + the antisymmetry permutation term  = F^k_j t2^{ab}_{ik}
    // addition A term: (params.pccsd_alpha - 1.0)/2.0 * A
    {
      const double A_coeff = (params.pccsd_alpha - 1.0)/2.0; //

      // allocate space for F
      dpdfile2 fIJ, FMI;
      dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
      dpd_file2_copy(&fIJ, CC_OEI, "FMI");
      dpd_file2_close(&fIJ);
      // compute f
      dpdbuf4 D, tIjAb, newtIjAb, Z;
      dpd_file2_init(&FMI, CC_OEI, 0, 0, 0, "FMI");
      dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
      dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
      dpd_contract442(&D, &tIjAb, &FMI, 0, 0, 1, 0.0);
      dpd_buf4_close(&tIjAb);
      dpd_buf4_close(&D);
      dpd_file2_close(&FMI);
      // contract f with t2
      dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
      dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
      dpd_file2_init(&FMI, CC_OEI, 0, 0, 0, "FMIt");
      dpd_contract244(&FMI, &tIjAb, &Z, 0, 0, 0, 1, 0);
      dpd_file2_close(&FMI);
      dpd_buf4_close(&tIjAb);
      dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
      dpd_buf4_axpy(&Z, &newtIjAb, -1.0 * A_coeff);
      dpd_buf4_close(&newtIjAb);
      dpd_buf4_sort_axpy(&Z, CC_TAMPS, qpsr, 0, 5, "New tIjAb", -1.0 * A_coeff);
      dpd_buf4_close(&Z);
    }

    //
    // B = 0.25 * v^{kl}^{cd} * t2^{cd}_{ij} t2^{ab}_{kl}  = 0.25 * w^{kl}_{ij} * t2^{ab}_{kl} (B is the term from the pCCSD paper)
    // the addition B term from pCCSD is: (pccsd_alpha -1) * B
    {
      const double B_coeff = params.pccsd_alpha - 1.0; //
      dpdbuf4 A, D, wMnIj, tIjAb, newtIjAb;

      // allocate space for w
      dpd_buf4_init(&A, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
      dpd_buf4_copy(&A, CC_HBAR, "wMnIj");
      dpd_buf4_close(&A);
      // compute w
      dpd_buf4_init(&wMnIj, CC_HBAR, 0, 0, 0, 0, 0, 0, "wMnIj");
      dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
      dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
      dpd_contract444(&D, &tIjAb, &wMnIj, 0, 0, 1.0, 0.0);
      dpd_buf4_close(&tIjAb);
      dpd_buf4_close(&D);
      dpd_buf4_close(&wMnIj);
      // contract w and t2 to produce B
      dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
      dpd_buf4_init(&wMnIj, CC_HBAR, 0, 0, 0, 0, 0, 0, "wMnIj");
      dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
      dpd_contract444(&wMnIj, &tIjAb, &newtIjAb, 1, 1, B_coeff, 1.0);
      dpd_buf4_close(&tIjAb);
      dpd_buf4_close(&wMnIj);
      dpd_buf4_close(&newtIjAb);
    }

  }

  if (params.pccsd_beta != 1.0) {

  }

}
}} // namespace psi::ccenergy
