#include <cstdlib>
#include <iostream>
#include <cmath>

#include <liboptions/liboptions.h>

#include "scf.h"

extern FILE* outfile;

namespace psi{ namespace mcscf{

void SCF::iterate_scf_equations()
{
  fprintf(outfile,"\n\n  ===========================================================================");
  fprintf(outfile,"\n         Cycle          Energy               D(Energy)            DIIS");
  fprintf(outfile,"\n  ===========================================================================");


  bool   converged  = false;
  int    cycle      = 0;
  double old_energy = 0.0;
  double new_energy = 0.0;
  while(!converged){
    C_T = C;
    C_T.transpose();

    // Compute the density matrix
    density_matrix();

    // Form F in the AO basis
    construct_F();

    // Compute energy
    old_energy = new_energy;
    new_energy = energy(cycle,old_energy);
    double delta_energy = new_energy - old_energy;

    // Transform F to the MO basis
    transform(Fc,Fc_t,C);
    if(reference == rohf)
      transform(Fo,Fo_t,C);
    if(reference == tcscf){
      for(int I = 0 ; I < nci; ++I)
        transform(Ftc[I],Ftc_t[I],C);
      transform(Favg,Favg_t,C);
    }

    construct_Feff(cycle);

    if(use_diis){
      diis(cycle);
      if(cycle >= ndiis){
        Feff_oAO.diagonalize(C_t,epsilon);
        C.multiply(false,false,S_sqrt_inv,C_t);
        C_T = C;
        C_T.transpose();
      }else{
        Feff_t.diagonalize(C_t,epsilon);
        T.multiply(false,false,C,C_t);
        C = T;
      }
    }else{
      Feff_t.diagonalize(C_t,epsilon);
      T.multiply(false,false,C,C_t);
      C = T;
    }


    if( fabs(log10(fabs(delta_energy))) > options_get_int("CONVERGENCE") ){
      if(reference == tcscf){
        if(2.0 * fabs(log10(norm_ci_grad)) > options_get_int("CONVERGENCE") )
          converged = true;
      }else{
        converged = true;
      }
    }

    if(cycle>options_get_int("MAX_ITERATIONS")){
      fprintf(outfile,"\n\n  The calculation did not converge in %d cycles",options_get_int("MAX_ITERATIONS"));
      fprintf(outfile,"\n  Quitting PSIMRCC\n");
      fflush(outfile);
      exit(1);
    }

    cycle++;
  }

  fprintf(outfile,"\n  ===========================================================================");
  fprintf(outfile,"\n\n  @SCF@     E(SCF) =  %20.12f",new_energy);
  
  fprintf(outfile,"\n\n  End of SCF");
  fflush(outfile);
}

}} /* End Namespaces */
