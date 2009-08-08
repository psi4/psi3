/**
 *  @file ccmrcc_pert_triples.cpp
 *  @ingroup (PSIMRCC)
 *  @brief Computes the (T) correction
*/

#include <libmoinfo/libmoinfo.h>
#include <libutil/libutil.h>
#include <libchkpt/chkpt.hpp>

#include "mrcc.h"
#include "mrccsd_t.h"


extern FILE* outfile;

namespace psi{ namespace psimrcc{

void CCMRCC::compute_perturbative_triples()
{
  Timer timer;
  fprintf(outfile,"\n\n\n  Computing (T) correction");
  fflush(outfile);
  
  h_eff.set_eigenvalue(current_energy);
  h_eff.set_matrix(Heff,moinfo->get_nrefs());
  h_eff.set_right_eigenvector(right_eigenvector,moinfo->get_nrefs());
  h_eff.set_left_eigenvector(left_eigenvector,moinfo->get_nrefs());

  MRCCSD_T mrccsd_t(&h_eff);

  current_energy = h_eff.expectation_value();
  _default_chkpt_lib_->wt_etot(current_energy);

  fprintf(outfile,"\n\n%6c* Mk-MRCCSD(T) total energy   =    %20.12f",' ',current_energy);
  fprintf(outfile,"\n\n  Timing for triples:             %20.6f s",timer.get());
  fflush(outfile);
//
//  for(int mu = 0; mu < moinfo->get_ref_size(UniqueRefs); ++mu){
//    int unique_mu = moinfo->get_ref_number(mu,UniqueRefs);
//    E4_all[unique_mu] = E4[mu];
//  }
//
//  for(int mu = 0; mu < moinfo->get_ref_size(AllRefs); ++mu){
//    int unique_mu = moinfo->get_ref_number(mu);
//    Heff[mu][mu] += E4_all[unique_mu];
//  }
}


}}  /* End Namespaces */
