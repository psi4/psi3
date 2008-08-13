#include <cstdlib>

#include <libmoinfo/libmoinfo.h>
#include <libutil/libutil.h>
#include <liboptions/liboptions.h>
#include <libchkpt/chkpt.h>

#include "blas.h"
#include "mp2_ccsd.h"
#include "debugging.h"
#include "matrix.h"
#include "sort.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{

using namespace std;

MP2_CCSD::MP2_CCSD():CCManyBody()
{
  triples_type = pt2;
  add_matrices();
}

MP2_CCSD::~MP2_CCSD()
{
}

void MP2_CCSD::compute_mp2_ccsd_energy()
{
  print_method("  MP2_CCSD\n    Using the DPD Library");

  generate_integrals();
  generate_denominators();
  compute_reference_energy();

  build_offdiagonal_F();

  blas->diis_add("t2[oO][vV]{u}","t2_delta[oO][vV]{u}");

  // Start the MP2 cycle
  bool converged = false;
  int  cycle     = 0;
  delta_energy = 0.0;
  while(!converged){
    build_mp2_t2_iJaB_amplitudes();
    blas->diis_save_t_amps(cycle);
    blas->diis(cycle,delta_energy,DiisEachCycle);

    blas->solve("t2[oo][vv]{u}  = t2[oO][vV]{u}");
    blas->solve("t2[oo][vv]{u} += #2134# - t2[oO][vV]{u}");
    blas->solve("t2[OO][VV]{u}  = t2[oo][vv]{u}");

    synchronize_amps(); // TODO: make this more efficient
    build_tau();  

    current_energy = compute_energy();
    delta_energy   = current_energy - old_energy;
    old_energy = current_energy;

    fprintf(outfile,"\n    @MP2      %5d   %20.15f  %11.4e",cycle,current_energy,delta_energy);

    if(fabs(log10(fabs(delta_energy))) > options_get_int("E_CONVERGENCE")){
      converged=true;
    }

    cycle++;
    fflush(outfile);
  }

  blas->diis_add("t1[o][v]{u}","t1_delta[o][v]{u}");

  // Start the MP2-CCSD cycle
  converged = false;
  cycle     = 0;
  delta_energy = 0.0;
  while(!converged){
    // These two go together before updating any other intermediate
    build_F_intermediates();
    build_W_intermediates();
    build_Z_intermediates();

    build_amplitudes();
    blas->diis_save_t_amps(cycle);
    blas->diis(cycle,delta_energy,DiisEachCycle);

    blas->solve("t2[oo][vv]{u}  = t2[oO][vV]{u}");
    blas->solve("t2[oo][vv]{u} += #2134# - t2[oO][vV]{u}");
    blas->solve("t2[OO][VV]{u}  = t2[oo][vv]{u}");
    blas->solve("t1[O][V]{u} = t1[o][v]{u}");

    synchronize_amps();
    build_tau();  

    current_energy = compute_energy();

    delta_energy = current_energy-old_energy;
    if(fabs(log10(fabs(delta_energy))) > options_get_int("E_CONVERGENCE")){
      converged=true;
    }
    fprintf(outfile,"\n    @MP2-CCSD %5d   %20.15f  %11.4e",cycle,current_energy,delta_energy);
    old_energy=current_energy;

    if(cycle>options_get_int("MAX_ITERATIONS")){
      fprintf(outfile,"\n\n\tThe calculation did not converge in %d cycles\n\tQuitting PSIMRCC\n",options_get_int("MAX_ITERATIONS"));
      fflush(outfile);
      exit(1);
    }
    cycle++;
    fflush(outfile);
  }
  fprintf(outfile,"\n  ------------------------------------------------------------------------------");
  
  fprintf(outfile,"\n\n    @MP2-CCSD@  =%25.15f\n",current_energy);

  fflush(outfile);
}

double MP2_CCSD::compute_energy()
{
  // Compute the energy using a simple UCCSD energy expression
  blas->solve("Eaa{u}   = t1[o][v]{u} . fock[o][v]{u}");
  blas->solve("Ebb{u}   = t1[O][V]{u} . fock[O][V]{u}");

  blas->solve("Eaaaa{u} = 1/4 tau[oo][vv]{u} . <[oo]:[vv]>");
  blas->solve("Eabab{u} =     tau[oO][vV]{u} . <[oo]|[vv]>");
  blas->solve("Ebbbb{u} = 1/4 tau[OO][VV]{u} . <[oo]:[vv]>");

  blas->solve("EPT2{u}  = Eaa{u} + Ebb{u} + Eaaaa{u} + Eabab{u} + Ebbbb{u} + ERef{u}");

  return(blas->get_scalar("EPT2",0));
}

void MP2_CCSD::read_mp2_ccsd_integrals()
{
  START_TIMER(1,"Reading the integrals required by MP2-CCSD");

  // CCSort reads the one and two electron integrals
  // and creates the Fock matrices
  sorter = new CCSort(out_of_core_sort);

  END_TIMER(1);
}

void MP2_CCSD::build_amplitudes()
{
  // These are required by the t1 amplitude equations
  build_t1_ia_amplitudes();
  build_t1_IA_amplitudes();
  build_t2_iJaB_amplitudes();
  build_t2_ijab_amplitudes();
  build_t2_IJAB_amplitudes();
}

}} /* End Namespaces */
