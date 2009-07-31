#include <cmath>
#include <cstdio>
#include <libchkpt/chkpt.h>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libipv1/ip_lib.h>
#include <physconst.h>
#include <libpsio/psio.h>
#include <psifiles.h>

#define EXTERN
#include "opt.h"
#undef EXTERN
#include "cartesians.h"
#include "internals.h"
#include "salc.h"
#include "bond_lengths.h"

namespace psi { namespace optking {

// angle bends cannot go < 0

void check_zero_angles(internals &simples, salc_set &symm, double *dq) {
  int i, j, dim, simple, intco_type, sub_index, sub_index2;
  double coeff, prefactor, dq_simple, val;
  
  dim = symm.get_num();

  // find bend angles with range {0,180}
  for (i=0; i<dim; ++i) { // loop over symm SALCS vectors
    prefactor = symm.get_prefactor(i);
    for (j=0;j<symm.get_length(i);++j) { // loop over simples in each salc
      coeff = symm.get_coeff(i,j);
      dq_simple = prefactor * coeff * dq[i]; // change in value of primitive

      // dq is in radians
      simple = symm.get_simple(i,j);
      simples.locate_id(simple,&intco_type,&sub_index,&sub_index2);

      if (intco_type == BEND_TYPE) {
        val = simples.bend.get_val(sub_index) / 180 * _pi;
        if (val + dq_simple < 0.0) {
          punt("Bond angle passing through zero. Try different torsions or higher symmetry.\n");
          exit(PSI_RETURN_FAILURE);
        }
      }
      else if ((intco_type == FRAG_TYPE) && (sub_index2 == 1)) {
        val = simples.frag.get_val(sub_index,sub_index2) / 180 * _pi;
        if (val + dq_simple < 0.0) {
          punt("Theta-A angle passing through zero. Try adding 180 to tau and chi-A or higher symmetry.\n");
          exit(PSI_RETURN_FAILURE);
        }
      }
      else if ((intco_type == FRAG_TYPE) && (sub_index2 == 2)) {
        val = simples.frag.get_val(sub_index,sub_index2) / 180 * _pi;
        if (val + dq_simple < 0.0) {
          punt("Theta-B angle passing through zero. Try adding 180 to tau and chi-B or higher symmetry.\n");
          exit(PSI_RETURN_FAILURE);
        }
      }
    }
  }
  return;
}

}}
