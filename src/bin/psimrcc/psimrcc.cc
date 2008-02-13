/***************************************************************************
 *  PSIMRCC
 *  Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/

/**
 *  @defgroup PSIMRCC psimrcc: a code for SR/MRCC computations
 *  @file psimrcc.cpp
 *  @ingroup PSIMRCC
 *  @brief Contains main() and global variables
*/

// Version 0.2   is an attempt to introduce matrix-matrix multuplication
// Version 0.2.1 introduces faster reindexing
//               CCPair was turned into CCIndex
//               fixed some bugs in the Fock matrix routines (make_fock_two)
//               fixed some bugs in the F intermediates routines (F_AB,F_mi,F_MI,F'_mi,F'_MI)
// Version 0.2.2 introduces restricted 4 indices virtual quantities
//               the "check" feature now allows to determine the memory required for a computation
//               parallel reindexing improved
//               a new out of core algorithm for the contraction with integrals
//               a primitive management of the in-core and out-of-core quantities
// Version 0.2.3 avoid the computation of open shell determinants that are related by spin flip
//               spin-adapted CCSD equations
// Version 0.5   Compute MkPT2 energies
//               The model space is manipulated by MOInfo
// Version 0.6   The CCBLAS class was completely rewritten (still debugging)
//               A partial test suite is being compiled
//               - The CCMatrix class was modified to facilitate
//                 in-core/out-of-core computations
// Version 0.6.2 Coded the MkCCSDT-1a equations and tested for closed- and open-shells
// Version 0.6.3 Out-of-core DIIS
// Version 0.6.4 Implements proper Mk-MRPT2


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "calculation_options.h"
#include "blas.h"
#include "sort.h"
#include "mrcc.h"
#include "idmrpt2.h"
#include "transform.h"
#include "debugging.h"
#include "moinfo.h"
#include "psimrcc.h"
#include "utilities.h"

namespace psi{ namespace psimrcc{

using namespace std;

void run_psimrcc()
{
  blas   = new CCBLAS();
  trans  = new CCTransform();

  if(options->get_str_option("CORR_WFN")=="PT2"){
    mrpt2();
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

  if(options->get_str_option("CORR_ANSATZ")=="SR")
    mrcc.compute_ccsd_energy();
  if(options->get_str_option("CORR_ANSATZ")=="MK")
    mrcc.compute_mkccsd_energy();
  if(options->get_str_option("CORR_ANSATZ")=="BW")
      mrcc.compute_bwccsd_energy();
  if(options->get_str_option("CORR_ANSATZ")=="APBW")
      mrcc.compute_apbwccsd_energy();
}

/*!
 * Runs a MRPT2 computation
 */
void mrpt2()
{
  // Initialize the mp2 module (integrals,fock matrix(ces),denominators)
  IDMRPT2        idmrpt2;

  // Compute the initial amplitudes and MP2 energy
  idmrpt2.compute_mrpt2_energy();

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
