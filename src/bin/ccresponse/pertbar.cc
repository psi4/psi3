/*! \file
    \ingroup CCRESPONSE
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstring>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

/* pertbar(): Constructs the components of the similarity-transformed
** one-electron property operator. 
**
** NB: For antisymmetric quantities, we must be careful in
** contractions involving the bare MO integrals.  This only affects
** the leading term of XXBAR_IA (which, I admit, is poorly named and
** probably should really be called XXBAR_AI), where I copy the XXIA
** integrals into place with a -1 sign.
**
** TDC, 10/05
*/

void pertbar(char *pert, int irrep_x, int irrep_y, int irrep_z, int anti)
{
  int irrep, i;
  dpdfile2 fbar1, f, t1, z;
  dpdbuf4 t2, fbar2;
  char lbl[32], prefix1[32], prefix2[32], prefix3[32];

  for(i=0; i < 3; i++) {
    if(i==0) {
      irrep = irrep_x;
      sprintf(prefix1, "%s_%1s", pert, "X");
      sprintf(prefix2, "%sBAR_%1s", pert, "X");
      sprintf(prefix3, "z_%1s", "X");
    }
    else if(i==1) {
      irrep = irrep_y;
      sprintf(prefix1, "%s_%1s", pert, "Y");
      sprintf(prefix2, "%sBAR_%1s", pert, "Y");
      sprintf(prefix3, "z_%1s", "Y");
    }
    else if(i==2) {
      irrep = irrep_z;
      sprintf(prefix1, "%s_%1s", pert, "Z");
      sprintf(prefix2, "%sBAR_%1s", pert, "Z");
      sprintf(prefix3, "z_%1s", "Z");
    }

    /** XXBAR_ME **/
    sprintf(lbl, "%s_IA", prefix1);
    dpd_file2_init(&f, CC_OEI, irrep, 0, 1, lbl);
    sprintf(lbl, "%s_ME", prefix2);
    dpd_file2_copy(&f, CC_OEI, lbl);
    dpd_file2_close(&f);

    /** XXBAR_MI **/
    sprintf(lbl, "%s_IJ", prefix1);
    dpd_file2_init(&f, CC_OEI, irrep, 0, 0, lbl);
    sprintf(lbl, "%s_MI", prefix2);
    dpd_file2_copy(&f, CC_OEI, lbl);
    dpd_file2_close(&f);

    sprintf(lbl, "%s_MI", prefix2);
    dpd_file2_init(&fbar1, CC_OEI, irrep, 0, 0, lbl);
    sprintf(lbl, "%s_IA", prefix1);
    dpd_file2_init(&f, CC_OEI, irrep, 0, 1, lbl);
    dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&f, &t1, &fbar1, 0, 0, 1, 1);
    dpd_file2_close(&t1);
    dpd_file2_close(&f);
    dpd_file2_close(&fbar1);

    /** XXBAR_AE **/
    sprintf(lbl, "%s_AB", prefix1);
    dpd_file2_init(&f, CC_OEI, irrep, 1, 1, lbl);
    sprintf(lbl, "%s_AE", prefix2);
    dpd_file2_copy(&f, CC_OEI, lbl);
    dpd_file2_close(&f);

    sprintf(lbl, "%s_AE", prefix2);
    dpd_file2_init(&fbar1, CC_OEI, irrep, 1, 1, lbl);
    sprintf(lbl, "%s_IA", prefix1);
    dpd_file2_init(&f, CC_OEI, irrep, 0, 1, lbl);
    dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&t1, &f, &fbar1, 1, 1, -1, 1);
    dpd_file2_close(&t1);
    dpd_file2_close(&f);
    dpd_file2_close(&fbar1);

    /** XXBAR_IA **/
    sprintf(lbl, "%s_IA", prefix1);
    dpd_file2_init(&f, CC_OEI, irrep, 0, 1, lbl);
    sprintf(lbl, "%s_IA", prefix2);
    dpd_file2_copy(&f, CC_OEI, lbl);
    dpd_file2_close(&f);

    sprintf(lbl, "%s_IA", prefix2);
    dpd_file2_init(&fbar1, CC_OEI, irrep, 0, 1, lbl);
    if(anti) dpd_file2_scm(&fbar1, -1);

    sprintf(lbl, "%s_AB", prefix1);
    dpd_file2_init(&f, CC_OEI, irrep, 1, 1, lbl);
    dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&t1, &f, &fbar1, 0, 0, 1, 1);
    dpd_file2_close(&t1);
    dpd_file2_close(&f);

    sprintf(lbl, "%s_IJ", prefix1);
    dpd_file2_init(&f, CC_OEI, irrep, 0, 0, lbl);
    dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&f, &t1, &fbar1, 1, 1, -1, 1);
    dpd_file2_close(&t1);
    dpd_file2_close(&f);

    sprintf(lbl, "%s_IA", prefix1);
    dpd_file2_init(&f, CC_OEI, irrep, 0, 1, lbl);
    dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "2 tIAjb - tIBja");
    dpd_contract422(&t2, &f, &fbar1, 0, 0, 1, 1);
    dpd_buf4_close(&t2);
    dpd_file2_close(&f);

    sprintf(lbl, "%s_MI", prefix3);
    dpd_file2_init(&z, CC_TMP0, irrep, 0, 0, lbl);
    dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
    sprintf(lbl, "%s_IA", prefix1);
    dpd_file2_init(&f, CC_OEI, irrep, 0, 1, lbl);
    dpd_contract222(&f, &t1, &z, 0, 0, 1, 0);
    dpd_file2_close(&f);
    dpd_contract222(&z, &t1, &fbar1, 1, 1, -1, 1);
    dpd_file2_close(&t1);
    dpd_file2_close(&z);

    dpd_file2_close(&fbar1);

    /** LBAR_MbIj **/
    sprintf(lbl, "%s_MbIj", prefix2);
    dpd_buf4_init(&fbar2, CC_LR, irrep, 10, 0, 10, 0, 0, lbl);
    sprintf(lbl, "%s_IA", prefix1);
    dpd_file2_init(&f, CC_OEI, irrep, 0, 1, lbl);
    dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_contract244(&f, &t2, &fbar2, 1, 2, 0, 1, 0);
    dpd_buf4_close(&t2);
    dpd_file2_close(&f);
    dpd_buf4_close(&fbar2);

    /** LBAR_IjAb **/
    sprintf(lbl, "%s_IjAb", prefix2);
    dpd_buf4_init(&fbar2, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);

    sprintf(lbl, "%s_AB", prefix1);
    dpd_file2_init(&f, CC_OEI, irrep, 1, 1, lbl);
    dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_contract424(&t2, &f, &fbar2, 3, 1, 0, 1, 0);
    dpd_contract244(&f, &t2, &fbar2, 1, 2, 1, 1, 1);
    dpd_buf4_close(&t2);
    dpd_file2_close(&f);

    sprintf(lbl, "%s_IJ", prefix1);
    dpd_file2_init(&f, CC_OEI, irrep, 0, 0, lbl);
    dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_contract424(&t2, &f, &fbar2, 1, 0, 1, -1, 1);
    dpd_contract244(&f, &t2, &fbar2, 0, 0, 0, -1, 1);
    dpd_buf4_close(&t2);
    dpd_file2_close(&f);

    sprintf(lbl, "%s_AE", prefix3);
    dpd_file2_init(&z, CC_TMP0, irrep, 1, 1, lbl);
    dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
    sprintf(lbl, "%s_IA", prefix1);
    dpd_file2_init(&f, CC_OEI, irrep, 0, 1, lbl);
    dpd_contract222(&t1, &f, &z, 1, 1, -1, 0);
    dpd_file2_close(&f);
    dpd_file2_close(&t1);
    dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_contract424(&t2, &z, &fbar2, 3, 1, 0, 1, 1);
    dpd_contract244(&z, &t2, &fbar2, 1, 2, 1, 1, 1);
    dpd_buf4_close(&t2);
    dpd_file2_close(&z);

    sprintf(lbl, "%s_MI", prefix3);
    dpd_file2_init(&z, CC_TMP0, irrep, 0, 0, lbl); /* generated above */
    dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_contract424(&t2, &z, &fbar2, 1, 0, 1, -1, 1);
    dpd_contract244(&z, &t2, &fbar2, 0, 0, 0, -1, 1);
    dpd_buf4_close(&t2);
    dpd_file2_close(&z);

    dpd_buf4_close(&fbar2);
  }
}

}} // namespace psi::ccresponse
