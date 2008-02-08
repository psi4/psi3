/***************************************************************************
 *  PSIMRCC : Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/

#include "calculation_options.h"
#include "moinfo.h"
#include "mrcc.h"
#include "matrix.h"
#include "blas.h"
#include "debugging.h"
#include "utilities.h"
#include "algebra_interface.h"

namespace psi{ namespace psimrcc{

using namespace std;

CCMRCC::CCMRCC()
{
  triples_type = ccsd;
  triples_coupling_type = cubic;
  ap_correction   = false; // Set tu true when computing the a posteriori correction
  current_energy  =  0.0;
  old_energy      = 10.0;

  // Add the matrices that will store the intermediates
  add_matrices();

  // Generate the Fock matrices, Integrals and Denominators
  generate_integrals();
  generate_denominators();
  if(triples_type>ccsd)
    generate_triples_denominators();

  compute_reference_energy();

  DEBUGGING(1,
    blas->print_memory();
  )

  // Parse the CORR_WFN parameter

  vector<string> theory_levels = split("PT2 CCSD CCSD_T CCSDT-1A CCSDT-1B CCSDT-2 CCSDT-3 CCSDT");
  for(int i=0;i<theory_levels.size();++i){
    if(options->get_str_option("CORR_WFN")==theory_levels[i])
      triples_type = TriplesType(i);
  }

  // Parse the COUPLING parameter

  vector<string> coupling_levels = split("NONE LINEAR QUADRATIC CUBIC");
  for(int i=0;i<coupling_levels.size();++i){
    if(options->get_str_option("COUPLING")==coupling_levels[i]){
      triples_coupling_type = TriplesCouplingType(i);
    }
  }
}

CCMRCC::~CCMRCC()
{
}

}} /* End Namespaces */