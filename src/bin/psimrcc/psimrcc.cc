/***************************************************************************
 *  PSIMRCC
 *  Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/

/**
 *  @defgroup PSIMRCC PSIMRCC is a code for SR/MRCC computations
 *  @file psimrcc.cpp
 *  @ingroup (PSIMRCC)
 *  @brief Contains main() and global variables
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <liboptions/liboptions.h>
#include "blas.h"
#include "sort.h"
#include "mp2_ccsd.h"
#include "mrcc.h"
#include "transform.h"
#include "debugging.h"
#include <libmoinfo/libmoinfo.h>
#include "psimrcc.h"
#include <libutil/libutil.h>

namespace psi{ namespace psimrcc{

using namespace std;

void run_psimrcc()
{
  blas   = new CCBLAS();
  trans  = new CCTransform();

  if(options_get_str("CORR_WFN")=="MP2-CCSD"){
    mp2_ccsd();
  }else{
    mrccsd();
  }

  delete sorter;
  delete trans;
  delete blas;
}

/*!
 * Runs a MRPT2 and a MRCCSD computation
 * @todo move this code in the CCMRCC class
 */
void mrccsd()
{
  // Initialize the mp2 module (integrals,fock matrix(ces),denominators)
  CCMRCC        mrcc;

  if(options_get_str("CORR_ANSATZ")=="SR")
    mrcc.compute_ccsd_energy();
  if(options_get_str("CORR_ANSATZ")=="MK")
    mrcc.compute_mkccsd_energy();
  if(options_get_str("CORR_ANSATZ")=="BW")
      mrcc.compute_bwccsd_energy();
  if(options_get_str("CORR_ANSATZ")=="APBW")
      mrcc.compute_apbwccsd_energy();
}

/*!
 * Runs a CCSD_MP2 computation
 */
void mp2_ccsd()
{
  // Initialize the mp2 module (integrals,fock matrix(ces),denominators)
  MP2_CCSD        mp2_ccsd;

  // Compute the initial amplitudes and CCSD_MP2 energy
  mp2_ccsd.compute_mp2_ccsd_energy();

  DEBUGGING(1,
    blas->print_memory();
  )
}

/*!
 * Runs a integral transformation
 * @todo CCTransform is still unused in the code
 */
void transform_integrals()
{
//   CCTransform transf;
//   transf.read_so_integrals();
}

}} /* End Namespaces */
