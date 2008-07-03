/*! \file
    \ingroup OPTKING
    \brief This function constructs the B matrix for a set of salcs
*/

# include <cmath>
#include <cstdio>
#include <libchkpt/chkpt.h>
#include <cstdlib>
#include <cstring>
#include <physconst.h>
#include <libciomr/libciomr.h>

#define EXTERN
#include "opt.h"
#undef EXTERN
#include "cartesians.h"
#include "internals.h"
#include "salc.h"

namespace psi { namespace optking {

double **compute_B(internals &simples,salc_set &symm) {
  int i,j,k,a,b,c,d, simple, intco_type, sub_index;
  int J,K,atom,xyz,cnt_frag=-1,nrefA,nrefB;
  double **B, coeff, prefactor, weight;

  B = block_matrix(symm.get_num(),3*optinfo.natom);

  for (i=0;i<symm.get_num();++i) {
    prefactor = symm.get_prefactor(i);
    for (j=0;j<symm.get_length(i);++j) {
      simple = symm.get_simple(i,j);
      coeff = symm.get_coeff(i,j);
      simples.locate_id(simple,&intco_type,&sub_index);
      if (intco_type == STRE_TYPE) {
        a = simples.stre.get_A(sub_index);
        b = simples.stre.get_B(sub_index);
        for (k=0;k<3;++k) {
          B[i][3*a+k] += prefactor * coeff * simples.stre.get_s_A(sub_index,k);
          B[i][3*b+k] += prefactor * coeff * simples.stre.get_s_B(sub_index,k);
        }
      }
      else if (intco_type == BEND_TYPE) {
        a = simples.bend.get_A(sub_index);
        b = simples.bend.get_B(sub_index);
        c = simples.bend.get_C(sub_index);
        for (k=0;k<3;++k) {
          B[i][3*a+k] += prefactor * coeff * simples.bend.get_s_A(sub_index,k);
          B[i][3*b+k] += prefactor * coeff * simples.bend.get_s_B(sub_index,k);
          B[i][3*c+k] += prefactor * coeff * simples.bend.get_s_C(sub_index,k);
        }
      }
      else if (intco_type == TORS_TYPE) {
        a = simples.tors.get_A(sub_index);
        b = simples.tors.get_B(sub_index);
        c = simples.tors.get_C(sub_index);
        d = simples.tors.get_D(sub_index);
        for (k=0;k<3;++k) {
          B[i][3*a+k] += prefactor * coeff * simples.tors.get_s_A(sub_index,k);
          B[i][3*b+k] += prefactor * coeff * simples.tors.get_s_B(sub_index,k);
          B[i][3*c+k] += prefactor * coeff * simples.tors.get_s_C(sub_index,k);
          B[i][3*d+k] += prefactor * coeff * simples.tors.get_s_D(sub_index,k);
        }
      }
      else if (intco_type == OUT_TYPE) {
        a = simples.out.get_A(sub_index);
        b = simples.out.get_B(sub_index);
        c = simples.out.get_C(sub_index);
        d = simples.out.get_D(sub_index);
        for (k=0;k<3;++k) {
          B[i][3*a+k] += prefactor * coeff * simples.out.get_s_A(sub_index,k);
          B[i][3*b+k] += prefactor * coeff * simples.out.get_s_B(sub_index,k);
          B[i][3*c+k] += prefactor * coeff * simples.out.get_s_C(sub_index,k);
          B[i][3*d+k] += prefactor * coeff * simples.out.get_s_D(sub_index,k);
        }
      }
      else if (intco_type == LIN_BEND_TYPE) {
        a = simples.lin_bend.get_A(sub_index);
        b = simples.lin_bend.get_B(sub_index);
        c = simples.lin_bend.get_C(sub_index);
        for (k=0;k<3;++k) {
          B[i][3*a+k] += prefactor * coeff * simples.lin_bend.get_s_A(sub_index,k);
          B[i][3*b+k] += prefactor * coeff * simples.lin_bend.get_s_B(sub_index,k);
          B[i][3*c+k] += prefactor * coeff * simples.lin_bend.get_s_C(sub_index,k);
        }
      }
      else if (intco_type == FRAG_TYPE) {
// for fragments with atoms < 3, fragments are assumed to appear in order and the
// interfragment distance coordinate must be included so that
// nref can guess the fragment number - a clean solution can await PSI4 :)
        if (simples.frag.get_J(sub_index) == 0) ++cnt_frag;
        if (optinfo.nref_per_fragment != NULL) {
          nrefA = optinfo.nref_per_fragment[cnt_frag];
          nrefB = optinfo.nref_per_fragment[cnt_frag+1];
        }
        else {
          nrefA = nrefB = 3; // which might be wrong - this data needs put into fragment class
        }
        for (a=0; a<simples.frag.get_A_natom(sub_index); ++a) { /* loop over ref atoms in A */
          atom   = simples.frag.get_A_atom(sub_index,a);        /* atom number of a'th atom in A */
          for (K=0;K<nrefA;++K) {  /* loop over reference atoms of A */
            weight = simples.frag.get_A_weight(sub_index,K,a);      /* weight of a'th atom in A for this K */
            for (xyz=0;xyz<3;++xyz)
              B[i][3*atom+xyz] += prefactor * coeff * weight * simples.frag.get_A_s(sub_index,3*K+xyz);
          }
        }
        for (b=0; b<simples.frag.get_B_natom(sub_index); ++b) { /* loop over ref atoms in B */
          atom   = simples.frag.get_B_atom(sub_index,b);        /* atom number of b'th atom in B */
          for (K=0;K<nrefB;++K) { /* loop over reference atoms of B*/
            weight = simples.frag.get_B_weight(sub_index,K,b);  /* weight of b'th atom in B for this K */
              for (xyz=0;xyz<3;++xyz)
                B[i][3*atom+xyz] += prefactor * coeff * weight * simples.frag.get_B_s(sub_index,3*K+xyz);
          }
        }
      }
    }
  }

  return B;
}

}} /* namespace psi::optking */
