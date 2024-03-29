#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <algorithm>
#include <cstdio>

#include <libmoinfo/libmoinfo.h>

#include "scf.h"

extern FILE* outfile;

using namespace std;

namespace psi{ namespace mcscf{

void SCF::print_eigenvectors_and_MO()
{
  typedef vector<std::pair<double, string> >::iterator vecstr_it;

  // Assumes the eigenvalues of some Fock operator
  // are in the SBlockVector epsilon
  vector<std::pair<double, string> > docc_evals;
  vector<std::pair<double, string> > actv_evals;
  vector<std::pair<double, string> > virt_evals;

  for(int h = 0; h < nirreps; ++h)
    for(int i = 0; i < docc[h]; ++i)
      docc_evals.push_back(
        make_pair(epsilon->get(h,i),moinfo_scf->get_irr_labs(h))
      );
  for(int h = 0; h < nirreps; ++h)
    for(int i = docc[h]; i < docc[h] + actv[h]; ++i)
      actv_evals.push_back(
        make_pair(epsilon->get(h,i),moinfo_scf->get_irr_labs(h))
      );
  for(int h = 0; h < nirreps; ++h)
    for(int i = docc[h] + actv[h]; i < sopi[h]; ++i)
      virt_evals.push_back(
        make_pair(epsilon->get(h,i),moinfo_scf->get_irr_labs(h))
      );

  sort(docc_evals.begin(),docc_evals.end());
  sort(actv_evals.begin(),actv_evals.end());
  sort(virt_evals.begin(),virt_evals.end());  


  fprintf(outfile,"\n\n  =========================================================================");
  fprintf(outfile,"\n  Eigenvalues (Eh)");
  fprintf(outfile,"\n  -------------------------------------------------------------------------");

  int print_nrows = 3;
  fprintf(outfile,"\n  Doubly occupied orbitals");
  int printed     = 0;
  int mo          = 1;
  for(vecstr_it it = docc_evals.begin(); it!= docc_evals.end(); ++it){
    fprintf(outfile,"%s  %5d %13.6f %3s",printed++ % print_nrows == 0 ? "\n" : "",
                                       mo++,
                                       it->first,
                                       it->second.c_str());
  }

  if(!actv_evals.empty()){
    fprintf(outfile,"\n  Active orbitals");
    printed     = 0;
    for(vecstr_it it = actv_evals.begin(); it!= actv_evals.end(); ++it){
      fprintf(outfile,"%s  %5d %13.6f %3s",printed++ % print_nrows == 0 ? "\n" : "",
                                        mo++,
                                        it->first,
                                        it->second.c_str());
    }
  }

  fprintf(outfile,"\n  Unoccupied orbitals");
  printed     = 0;
  for(vecstr_it it = virt_evals.begin(); it!= virt_evals.end(); ++it){
    fprintf(outfile,"%s  %5d %13.6f %3s",printed++ % print_nrows == 0 ? "\n" : "",
                                       mo++,
                                       it->first,
                                       it->second.c_str());
  }
  fprintf(outfile,"\n  =========================================================================\n");
}

}} // End namespace
