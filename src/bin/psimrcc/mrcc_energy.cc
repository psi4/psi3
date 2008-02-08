/***************************************************************************
 *  PSIMRCC : Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/
#include "moinfo.h"
#include "calculation_options.h"
#include "mrcc.h"
#include "matrix.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{

using namespace std;

void CCMRCC::print_mrccsd_energy(int cycle)
{
  double delta_energy = current_energy-old_energy;
  if(cycle==0){
    print_method("\tMultireference Coupled Cluster\n\t\tUsing the DPD Library");
    fprintf(outfile,"\n  ------------------------------------------------------------------------------");
    fprintf(outfile,"\n  @CC Cycle      Energy          Delta E    ||DeltaT1|| ||DeltaT2|| Timing  DIIS");
    fprintf(outfile,"\n  @CC           (Hartree)       (Hartree)                           (Sec)");
    fprintf(outfile,"\n  ------------------------------------------------------------------------------");
  }
  if(cycle>=0){
    fprintf(outfile,"\n  @CC %3d  %18.12f  %11.4e   %8.3e   %8.3e %7.0f",cycle,current_energy,delta_energy,delta_t1_amps,delta_t2_amps,total_time);

    if((fabs(log10(fabs(delta_energy))) > options->get_int_option("E_CONVERGENCE")) && (cycle!=0)){
      fprintf(outfile,"\n  ------------------------------------------------------------------------------");
      fprintf(outfile,"\n\n\t\t----------------------------------------------------------------");
        fprintf(outfile,"\n\t\t@CC@\t\t%s%s  Energy       =    %20.15f",options->get_str_option("CORR_ANSATZ").c_str(),options->get_str_option("CORR_WFN").c_str(),current_energy);
      fprintf(outfile,"\n\t\t----------------------------------------------------------------");
    }
  }else if(cycle==-1){
    fprintf(outfile,"\n\t-----------------------------------------------------------------------------------------");
    fprintf(outfile,"\n\n\t\t----------------------------------------------------------------");
      fprintf(outfile,"\n\t\t@CC2@\t\t%s%s  Energy     =    %20.15f",options->get_str_option("CORR_ANSATZ").c_str(),options->get_str_option("CORR_WFN").c_str(),current_energy);
    fprintf(outfile,"\n\t\t----------------------------------------------------------------");
    print_eigensystem(moinfo->get_nrefs(),Heff,eigenvector);
  }

  fflush(outfile);
}

}} /* End Namespaces */