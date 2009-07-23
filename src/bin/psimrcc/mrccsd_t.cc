#include "mrccsd_t.h"

namespace psi{ namespace psimrcc{

//using namespace std;

MRCCSD_T::MRCCSD_T(Hamiltonian* h_eff_) : h_eff(h_eff_)
{
  startup();
  compute();
}

MRCCSD_T::~MRCCSD_T()
{
  cleanup();
}

}} /* End Namespaces */
