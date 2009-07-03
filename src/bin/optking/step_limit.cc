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

// no internal coordinate value can change more than optinfo.step_limit
// dq is passed in in Angstroms or radians

void step_limit(internals &simples, salc_set &symm, double *dq) {

  int i, j, dim, max_i, simple, intco_type, sub_index, sub_index2;
  double scale = 1.0, coeff, prefactor, tval, max_dq = 0.0;
  
  dim = symm.get_num();

  for (i=0; i<dim; ++i) { // loop over symm vectors
    prefactor = symm.get_prefactor(i);
    for (j=0;j<symm.get_length(i);++j) { // loop through individual salc
      coeff = symm.get_coeff(i,j);
      tval = prefactor * coeff * dq[i];

      // determine if it is a strech, and we have to do bohr -> Angstrom conversion
      simple = symm.get_simple(i,j);
      simples.locate_id(simple,&intco_type,&sub_index,&sub_index2);

      if ( (intco_type == STRE_TYPE) || ((intco_type == FRAG_TYPE) && (sub_index2 == 5)) )
        tval /= _bohr2angstroms;

      if (fabs(tval) > max_dq) {
        max_dq = fabs(tval);
        max_i = i;
      }
    }
  }

  // avoid piddly symmetry breaking
  for (i=0;i<dim;++i)
    if (fabs(dq[i]) < MIN_DQ_STEP) dq[i] = 0.0;

  if (max_dq > optinfo.step_limit) {
    fprintf(outfile,"\nMaximum change in SALC %d exceeds STEP_LIMIT\n", max_i+1);
    scale = optinfo.step_limit / max_dq;
    fprintf(outfile,"Scaling displacements by %lf\n",scale);
    for (i=0;i<dim;++i)
      dq[i] *= scale;   
  }

  return;
}

}}
