#include <liboptions/liboptions.h>

#include "scf.h"

namespace psi{ namespace mcscf{

void SCF::construct_Feff(int cycle)
{
  Feff_t = Fc_t;

  if(options_get_bool("USE_FAVG")){
    if(cycle >= options_get_int("START_FAVG")){
      Feff_t = Favg_t;
    }
    if(cycle == options_get_int("START_FAVG")){
      fprintf(outfile,"\n  *** Switching from Fc to F_avg ***");
    }
  }

  // Set the diagonal blocks Fock 
//   for(int h =0; h < nirreps; ++h){
//     // Set the (virtual,virtual) blocks
//     for(int i = docc[h] + actv[h]; i < sopi[h]; ++i){
//       for(int j = docc[h] + actv[h]; j < sopi[h]; ++j){
//         
//         double element = 2.0 * Fo_t->get(h,i,j);
//         Feff_t->set(h,i,j,element);
//       }
//     }
//   }

  if(reference == rohf  && (cycle > turn_on_actv)){
    for(int h =0; h < nirreps; ++h){
      // Set the (closed,open) and (open,closed) blocks to 2(Fc-Fo)
      for(int i =0; i < docc[h]; ++i){
        for(int j = docc[h]; j < docc[h] + actv[h]; ++j){
          double element = 2.0 * ( Fc_t->get(h,i,j) - Fo_t->get(h,i,j) );
          Feff_t->set(h,i,j,element);
          Feff_t->set(h,j,i,element);
        }
      }
  
      // Set the (virtual,open) and (open,virtual) blocks to 2 Fo
      for(int i = docc[h] + actv[h]; i < sopi[h]; ++i){
        for(int j = docc[h]; j < docc[h] + actv[h]; ++j){
          double element = 2.0 * Fo_t->get(h,i,j);
          Feff_t->set(h,i,j,element);
          Feff_t->set(h,j,i,element);
        }
      }
    }
  }
  if(reference == tcscf && (cycle > turn_on_actv)){
    double lambda = 1.0;
    for(int I = 0 ; I < nci; ++I){
      int h = tcscf_sym[I];
      int i = tcscf_mos[I];
      // Set the (closed,tc) and (tc,closed) blocks to 2 (Fc - Ftc)
      for(int j = 0; j < docc[h]; ++j){
        double element = 2.0 * ( Fc_t->get(h,i,j) - Ftc_t[I]->get(h,i,j) );
        Feff_t->set(h,i,j,element);
        Feff_t->set(h,j,i,element);
      }
      // Set the (external,tc) and (tc,external) blocks to 2 Ftc
      for(int j = docc[h] + actv[h]; j < sopi[h]; ++j){
        double element = 2.0 * Ftc_t[I]->get(h,i,j);
        Feff_t->set(h,i,j,element);
        Feff_t->set(h,j,i,element);
      }
    }
    // Orbitals of the same symmetry
    for(int I = 0 ; I < nci; ++I){
      for(int J = I + 1 ; J < nci; ++J){
        if(tcscf_sym[I] == tcscf_sym[J]){
          int h = tcscf_sym[I];
          int i = tcscf_mos[I];
          int j = tcscf_mos[J];
          // Set the (tc,tc) and (tc,tc) blocks to 2 (Ftc_a - Ftc_b)
          double element = 2.0 * (Ftc_t[I]->get(h,i,j) - Ftc_t[J]->get(h,i,j));
          Feff_t->set(h,i,j,element);
          Feff_t->set(h,j,i,element);
        }
      }
    }
  }

}

}} /* End Namespaces */
