/***************************************************************************
 *  PSIMRCC : Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/

/**
 *  @file ccmrcc_compute.cpp
 *  @ingroup (PSIMRCC)
 *  @brief Contains all the methods to compute the energy
*/

#include <liboptions/liboptions.h>
#include "blas.h"
#include "mrcc.h"
#include "debugging.h"
#include <libmoinfo/libmoinfo.h>
#include <cstdlib>

namespace psi{ namespace psimrcc{

// typedef void(CCMRCC::*voidFunction)();

CCMRCC* CCMRCC::ptr = NULL;

void CCMRCC::compute_ccsd_energy()
{
  CCMRCC::ptr = this;
  compute_energy(update_amps_ccsd_wrapper);
}

void CCMRCC::compute_mkccsd_energy()
{
  CCMRCC::ptr = this;
  compute_energy(update_amps_mkccsd_wrapper);
}

void CCMRCC::compute_mkccsd_residual_energy()
{
  CCMRCC::ptr = this;
  compute_energy(update_amps_mkccsd_residual_wrapper);
}

void CCMRCC::compute_bwccsd_energy()
{
  CCMRCC::ptr = this;
  compute_energy(update_amps_bwccsd_wrapper);
}

void CCMRCC::compute_apbwccsd_energy()
{
  ap_correction = true;
  CCMRCC::ptr = this;
  compute_energy(update_amps_bwccsd_wrapper);
}

/**
 * This is a generic coupled cluster cycle. By specifying updater you can get all the falvours, single-reference, Mukherjee MRCC,...
 * @param (* updater)( ) the pointer to a CCMRCC static function that wrapps the updater function
 */
void CCMRCC::compute_energy(void(*updater)())
{
  blas->diis_add("t1[o][v]{u}","t1_delta[o][v]{u}");
  blas->diis_add("t1[O][V]{u}","t1_delta[O][V]{u}");
  blas->diis_add("t2[oo][vv]{u}","t2_delta[oo][vv]{u}");
  blas->diis_add("t2[oO][vV]{u}","t2_delta[oO][vV]{u}");
  blas->diis_add("t2[OO][VV]{u}","t2_delta[OO][VV]{u}");
  if(options_get_bool("DIIS_TRIPLES")){
    blas->diis_add("t3[ooo][vvv]{u}","t3_delta[ooo][vvv]{u}");
    blas->diis_add("t3[ooO][vvV]{u}","t3_delta[ooO][vvV]{u}");
    blas->diis_add("t3[oOO][vVV]{u}","t3_delta[oOO][vVV]{u}");
    blas->diis_add("t3[OOO][VVV]{u}","t3_delta[OOO][VVV]{u}");
  }

  Timer cc_timer;
  bool converged = false;
  // Start CC cycle
  int cycle = 0;
  while(!converged){
    diis_step = cycle % options_get_int("MAXDIIS");  

    zero_internal_amps();

    synchronize_amps();
    build_tau_intermediates();
    build_F_intermediates();
    build_W_intermediates();
    build_W_T3_intermediates();
    build_Z_intermediates();
    build_t1_amplitudes();
    build_t2_amplitudes();
    blas->compute();
    if(triples_type>ccsd_t)
      build_t1_amplitudes_triples();
    if(triples_type>ccsd_t)
      build_t2_amplitudes_triples();
    build_t3_amplitudes();

    converged=build_diagonalize_Heff(cycle,cc_timer.get());
    if(!converged){
      blas->diis_save_t_amps(cycle);
      updater();
      blas->diis(cycle,delta_energy,DiisCC);
    }

    if(cycle>options_get_int("MAX_ITERATIONS")){
      fprintf(outfile,"\n\n\tThe calculation did not converge in %d cycles\n\tQuitting PSIMRCC\n",options_get_int("MAX_ITERATIONS"));
      fflush(outfile);
      exit(1);
    }
    cycle++;
  }

  // Compute the apBWCCSD energy
  if(ap_correction){
    zero_internal_amps();

    synchronize_amps();

    build_tau_intermediates();
    build_F_intermediates();
    build_W_intermediates();
    build_Z_intermediates();

    build_t1_amplitudes();
    build_t2_amplitudes();

    update_amps();

    zero_internal_amps();

    synchronize_amps();

    build_tau_intermediates();
    build_F_intermediates();
    build_W_intermediates();
    build_Z_intermediates();

    build_t1_amplitudes();
    build_t2_amplitudes();

    zero_internal_amps();

    converged=build_diagonalize_Heff(-1,cc_timer.get());
  }

  DEBUGGING(1,
    blas->print_memory();
  );
  CCOperation::print_timing();
}

}} /* End Namespaces */
