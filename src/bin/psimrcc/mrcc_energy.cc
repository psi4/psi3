#include <libmoinfo/libmoinfo.h>
#include <liboptions/liboptions.h>
#include "mrcc.h"
#include "matrix.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{

using namespace std;

void CCMRCC::print_mrccsd_energy(int cycle)
{
  delta_energy = current_energy-old_energy;
  if(cycle==0){
    print_method("\tMultireference Coupled Cluster\n\t\tUsing the DPD Library");
    fprintf(outfile,"\n  ------------------------------------------------------------------------------");
    fprintf(outfile,"\n  @CC Cycle      Energy          Delta E    ||DeltaT1|| ||DeltaT2|| Timing  DIIS");
    fprintf(outfile,"\n  @CC           (Hartree)       (Hartree)                           (Sec)");
    fprintf(outfile,"\n  ------------------------------------------------------------------------------");
  }
  if(cycle>=0){
    fprintf(outfile,"\n  @CC %3d  %18.12f  %11.4e   %8.3e   %8.3e %7.0f",cycle,current_energy,delta_energy,delta_t1_amps,delta_t2_amps,total_time);

    
    if((fabs(log10(fabs(delta_energy))) > options_get_int("CONVERGENCE")) && (cycle!=0)){
      fprintf(outfile,"\n  ------------------------------------------------------------------------------");
      
      fprintf(outfile,"\n\n%6c* %s-MR%s total energy   =    %20.12f",' ',options_get_str("CORR_ANSATZ").c_str(),options_get_str("CORR_WFN").c_str(),current_energy);
//      fprintf(outfile,"\n\n%6c* Mk-MRCCSD total energy   = %20.12f\n",' ',current_energy);      
    }
  }else if(cycle==-1){
    fprintf(outfile,"\n\n%6c* Mk-MRCCSD total energy   = %20.12f\n",' ',current_energy);
    print_eigensystem(moinfo->get_nrefs(),Heff,eigenvector);
  }

  fflush(outfile);
}

}} /* End Namespaces */
