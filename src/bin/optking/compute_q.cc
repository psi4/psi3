/*! \file
    \ingroup OPTKING
    \brief Returns the values of salcs given the simple internals
    and salc_set the value of the simple internals must already be computed.
*/

#define EXTERN
#include "globals.h"
#undef EXTERN
#include "cartesians.h"
#include "simples.h"
#include "salc.h"

namespace psi { namespace optking {

double *compute_q(const simples_class & simples, const salc_set & symm) {
  int i, j, simple, sub_index, sub_index2;
  double *q, coeff, prefactor, tval;
  Intco_type intco_type;

  q = init_array(symm.get_num());

  /* q is build in angstroms and radians */
  for (i=0;i<symm.get_num();++i) {
    prefactor = symm.get_prefactor(i);
    for (j=0;j<symm.get_length(i);++j) {
      simple = symm.get_simple(i,j);
      coeff = symm.get_coeff(i,j);
      simples.locate_id(simple,&intco_type,&sub_index,&sub_index2);
      if (intco_type == STRE_TYPE) {
        q[i] += prefactor * coeff * simples.stre[sub_index].get_val();
      }
      else if (intco_type == BEND_TYPE) {
        q[i] += prefactor * coeff * simples.bend[sub_index].get_val()*_pi/180.0;
      }
      else if (intco_type == TORS_TYPE) {
        q[i] += prefactor * coeff * simples.tors[sub_index].get_val()*_pi/180.0;
      }
      else if (intco_type == OUT_TYPE) {
        q[i] += prefactor * coeff * simples.out[sub_index].get_val()*_pi/180.0;
      }
      else if (intco_type == LINB_TYPE) {
        tval = simples.linb[sub_index].get_val();
        q[i] += prefactor * coeff * tval*_pi/180.0;
      }
      else if (intco_type == FRAG_TYPE) {
        q[i] += prefactor * coeff * simples.frag[sub_index].get_val_A_or_rad(sub_index2);
      }
    }
  }

  return q;
}

}} /* namespace psi::optking */

