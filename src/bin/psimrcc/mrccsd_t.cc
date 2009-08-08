#include <liboptions/liboptions.h>

#include "mrccsd_t.h"

namespace psi{ namespace psimrcc{

MRCCSD_T::MRCCSD_T(Hamiltonian* h_eff_) : h_eff(h_eff_)
{
  startup();
  if(options_get_bool("RESTRICTED_TRIPLES"))
    compute_restricted();
  else
    compute();
}

MRCCSD_T::~MRCCSD_T()
{
  cleanup();
}

}} /* End Namespaces */
