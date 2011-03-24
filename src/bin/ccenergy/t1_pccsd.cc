/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"



namespace psi { namespace ccenergy {

/// so far only restricted version is implemented. This function is for the gamma factor term of pCCSD.
/// With respect to the document from Lee, the additional terms are (pccsd_gamma -1 ) * (E + F)
void t1_delta_pCCSD_build(void)
{
  int h,a,e,nirreps;
  int ma,fe,ef,m,f,M,A,Gm,Ga,Ge,Gf,Gma,nrows,ncols;
  double *X;
  dpdfile2 tIA, tia, newtIA;
  dpdfile2 FME, Fme;
  dpdfile2 fAB, fab, fIA, fia;
  dpdfile2 FAE, Fae;
  dpdfile2 FAEt, Faet;
  dpdbuf4 F_anti, F, D_anti, D;
  dpdbuf4 tIjAb;


  nirreps = moinfo.nirreps;

  const double E_fac = params.pccsd_gamma -1;
  const double F_fac = params.pccsd_gamma -1;



  if(params.ref == 0) { /** RHF **/
    /// for E term, first g is contracted to T2, then contracted to T1.
    /// the first contraction, tweaking of Fae_build
    dpd_file2_init(&newtIA, CC_OEI, 0, 0, 1, "New tIA");
    dpd_file2_init(&FAE, CC_OEI, 0, 1, 1, "FAE");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_contract442(&tIjAb, &D, &FAE, 3, 3, -1 * E_fac, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tIjAb);
    // then contract to T1
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&tIA, &FAE, &newtIA, 0, 0, 1, 1);
    dpd_file2_close(&FAE);
    dpd_file2_close(&tIA);


    /// for F term, first g is contracted to T1 (part of FME), then contracted to T2


    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");

    dpd_buf4_init(&D_anti, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij||ab>");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");

    dpd_dot13(&tIA, &D_anti, &FME, 0, 0, F_fac, 0.0);
    dpd_dot13(&tIA, &D, &FME, 0, 0, F_fac, 1.0);

    dpd_file2_close(&tIA);
    dpd_buf4_close(&D_anti);
    dpd_buf4_close(&D);

    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 10, 10, 10, 10, 0, "2 tIAjb - tIBja");
    dpd_contract422(&tIjAb, &FME, &newtIA, 0, 0, 1, 1);
    dpd_buf4_close(&tIjAb);
    dpd_file2_close(&FME);
    dpd_file2_close(&newtIA);
  }
}
}} // namespace psi::ccenergy






