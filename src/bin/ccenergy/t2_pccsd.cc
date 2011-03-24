
/*! \file
    \ingroup CCENERGY
    \brief implementation the addition terms needed for doubles residual of pCCSD of Nooijen;
*/
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include "Params.h"
#define EXTERN
#include "globals.h"
#include <cstring>
#include <assert.h>
#include "MOInfo.h"


namespace psi { namespace ccenergy {
  void WmbejT2(void);

void t2_delta_pCCSD_build(void)
{
  // only RHF implemented
  assert( params.ref == 0 );

  dpdbuf4 A_anti, A, Eijka, Eijka_anti, Eaijk, Eaijk_anti, D, D_anti;
  dpdbuf4 WMnIj, W;
  dpdbuf4 tIjAb, t2;
  dpdbuf4 Z, newtIJAB, newtijab, newtIjAb;
  dpdfile2 tIA, tia, FMIt, FAEt;


  if (params.pccsd_alpha != 1.0)
  {
    // A = -1/2 * v^{kl}_{cd}*  t2^{cd}_{jl} * t2^{ab}_{ik} + the antisymmetry permutation term  = F^k_j t2^{ab}_{ik} (refer to Nooijen pCCSD JCP paper)
    // additional A term: (params.pccsd_alpha - 1.0)/2.0 * A
    // the first contraction is obtained from modifying Fmi.cc by removing irrelevant terms (FMIt is reset);
    // the second contraction is obtained from modifying FmitT2.cc

    const double A_coeff = (params.pccsd_alpha - 1.0)/2.0;

    dpd_file2_init(&FMIt, CC_OEI, 0, 0, 0, "FMIt");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_contract442(&D, &tIjAb, &FMIt, 0, 0, 1.0, 0.0); // FMIt is reset by setting beta = 0.0
    dpd_buf4_close(&tIjAb);
    dpd_buf4_close(&D);

    dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    //dpd_file2_init(&FMIt, CC_OEI, 0, 0, 0, "FMIt");
    dpd_contract244(&FMIt, &tIjAb, &Z, 0, 0, 0, A_coeff, 0.0); // for easy implementation, we put A_coeff here;
    dpd_file2_close(&FMIt);
    dpd_buf4_close(&tIjAb);
    dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    dpd_buf4_axpy(&Z, &newtIjAb, -1);
    dpd_buf4_close(&newtIjAb);
    dpd_buf4_sort_axpy(&Z, CC_TAMPS, qpsr, 0, 5, "New tIjAb", -1);
    dpd_buf4_close(&Z);


    // B = 0.25 * v^{kl}_{cd} * t2^{cd}_{ij} t2^{ab}_{kl}  = 0.25 * w^{kl}_{ij} * t2^{ab}_{kl} (for definition of B, refer to Nooijen pCCSD JCP paper)
    // the additional B term from pCCSD is: (pccsd_alpha -1) * B
    // the first contraction is obtained from modifying Wmnij.cc by removing irrelevant terms (WMnIj is reset);
    // the second step contraction is obtained from modification of WmnijT2.cc

    const double B_coeff = params.pccsd_alpha - 1.0;

    dpd_buf4_init(&WMnIj, CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_contract444(&D, &tIjAb, &WMnIj, 0, 0, 1.0, 0.0); // Z = a X*Y + b Z; b set to 0 here, so WMnIj is reset
    dpd_buf4_close(&tIjAb);
    dpd_buf4_close(&D);
    dpd_buf4_close(&WMnIj);


    dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    dpd_buf4_init(&WMnIj, CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_contract444(&WMnIj, &tIjAb, &newtIjAb, 1, 1, B_coeff, 1.0); // B_coeff
    dpd_buf4_close(&tIjAb);
    dpd_buf4_close(&WMnIj);
    dpd_buf4_close(&newtIjAb);
  } //end of pccsd_alpha section


  if (params.pccsd_beta != 1.0)
  {

    // C = -0.5 * v^{kl}_{cd} * t2^{bd}_{kl} t2^{ac}_{ij} + ..  = -0.5 * w^b_c * t2^{ac}_{ij}
    // the additional C term from pCCSD is: (pccsd_beta -1) * C
    // the first contraction is obtained from modifying Fae.cc by removing irrelevant terms (FAEt is reset);
    // the second step contraction is obtained from modification of FaetT2.cc

      const double C_coeff = params.pccsd_beta -1;
      const double D_coeff = C_coeff;


      dpd_file2_init(&FAEt, CC_OEI, 0, 1, 1, "FAEt");
      dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
      dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
      dpd_contract442(&tIjAb, &D, &FAEt, 3, 3, -1, 0.0);  // reset FAEt
      dpd_buf4_close(&D);
      dpd_buf4_close(&tIjAb);
      dpd_file2_close(&FAEt);


      dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
      dpd_file2_init(&FAEt, CC_OEI, 0, 1, 1, "FAEt");
      dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "Zijab");
      dpd_contract424(&tIjAb, &FAEt, &Z, 3, 1, 0, C_coeff, 0); // C_coeff
      dpd_file2_close(&FAEt);
      dpd_buf4_close(&tIjAb);
      dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
      dpd_buf4_axpy(&Z, &newtIjAb, 1);
      dpd_buf4_close(&newtIjAb);
      dpd_buf4_sort_axpy(&Z, CC_TAMPS, qpsr, 0, 5, "New tIjAb", 1);
      dpd_buf4_close(&Z);

      // D = v^{kl}_{cd} * t2^{ac}_{ik}* t2^{bd}_{jl} + ..  = w^{al}_{id} * t2^{bd}_{jl}
      // the addition D term from pCCSD is: (pccsd_beta -1) * D
      // the first contraction is obtained from modifying Wmbej.cc by removing irrelevant terms (WMbeJ is reset);


        /*** ABAB ***/

        dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
        dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
        dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)");
        dpd_contract444(&D, &t2, &W, 0, 1, 0.5 * D_coeff, 0.0); // reset WMbEj
        dpd_buf4_close(&D);
        dpd_buf4_close(&t2);
        dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tjAIb");
        dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
        dpd_contract444(&D, &t2, &W, 0, 1, -0.5 * D_coeff, 1);
        dpd_buf4_close(&D);
        dpd_buf4_close(&t2);
        dpd_buf4_close(&W);


        /*** ABBA ***/


        //weird
        dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
        dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
        dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
        dpd_contract444(&D, &t2, &W, 0, 1, 0.5 * D_coeff, 0.0); // reset WMbeJ
        dpd_buf4_close(&D);
        dpd_buf4_close(&t2);
        dpd_buf4_close(&W);
        WmbejT2();

  }

}}}




