// returns the values of salcs given the simple internals and salc_set
// the value of the simple internals must already be computed

#if HAVE_CMATH
# include <cmath>
#else
# include <math.h>
#endif

extern "C" {
#include <stdio.h>
#include <libchkpt/chkpt.h>
#include <stdlib.h>
#include <string.h>
#include <libciomr/libciomr.h>
#include <physconst.h>
}

#define EXTERN
#include "opt.h"
#undef EXTERN
#include "cartesians.h"
#include "internals.h"
#include "salc.h"

double *compute_q(internals &simples,salc_set &symm) {
  int i, j, simple, intco_type, sub_index;
  double *q, coeff, prefactor, tval;

  q = init_array(symm.get_num());

  for (i=0;i<symm.get_num();++i) {
    prefactor = symm.get_prefactor(i);
    for (j=0;j<symm.get_length(i);++j) {
      simple = symm.get_simple(i,j);
      coeff = symm.get_coeff(i,j);
      simples.locate_id(simple,&intco_type,&sub_index);
      if (intco_type == STRE_TYPE) {
        q[i] += prefactor * coeff * simples.stre.get_val(sub_index);
      }
      else if (intco_type == BEND_TYPE) {
        q[i] += prefactor * coeff * simples.bend.get_val(sub_index)*_pi/180.0;
      }
      else if (intco_type == TORS_TYPE) {
        q[i] += prefactor * coeff * simples.tors.get_val(sub_index)*_pi/180.0;
      }
      else if (intco_type == OUT_TYPE) {
        q[i] += prefactor * coeff * simples.out.get_val(sub_index)*_pi/180.0;
      }
      else if (intco_type == LIN_BEND_TYPE) {
        tval = simples.lin_bend.get_val(sub_index);
        // if ( fabs(180.0-tval) < 0.001 ) tval = 180.0;
        q[i] += prefactor * coeff * tval*_pi/180.0;
      }
    }
  }

  return q;
}
