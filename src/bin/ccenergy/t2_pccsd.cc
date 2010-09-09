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
    dpdbuf4 A, D, wMnIj, tIjAb, newtIjAb;

    //
    // B = v_{}^{} t2_^ t2_^ = w_^ t2_^
    //
    // allocate space for w
    dpd_buf4_init(&A, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
    dpd_buf4_copy(&A, CC_HBAR, "wMnIj");
    dpd_buf4_close(&A);
    // compute w
    dpd_buf4_init(&wMnIj, CC_HBAR, 0, 0, 0, 0, 0, 0, "wMnIj");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_contract444(&D, &tIjAb, &wMnIj, 0, 0, 1, 1);
    dpd_buf4_close(&tIjAb);
    dpd_buf4_close(&D);
    dpd_buf4_close(&wMnIj);
    // contract w and t2 to produce B
    dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    dpd_buf4_init(&wMnIj, CC_HBAR, 0, 0, 0, 0, 0, 0, "wMnIj");
    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_contract444(&wMnIj, &tIjAb, &newtIjAb, 1, 1, (params.pccsd_alpha-1.0), 1.0);
    dpd_buf4_close(&tIjAb);
    dpd_buf4_close(&wMnIj);
    dpd_buf4_close(&newtIjAb);

  }

  if (params.pccsd_beta != 1.0) {

  }

}
}} // namespace psi::ccenergy
