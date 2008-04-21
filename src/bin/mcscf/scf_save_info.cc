#include <libchkpt/chkpt.h>

#include "scf.h"
#include "memory_manager.h"

namespace psi{ namespace mcscf{

void SCF::save_info()
{
  // Writes out the total number of irreps in the point group in which the molecule is being considered which have non-zero number of basis functions.
  int n_so_typs = 0;
  for(int h = 0; h < nirreps; ++h){
    if( docc[h] + actv[h] > 0 ) n_so_typs++;
  }
  chkpt_wt_nsymhf(n_so_typs);

  // Writes out the total number of molecular orbitals.
  chkpt_wt_nmo(nso); // TODO: find nmo

  // Writes out the dimensionality of ALPHA and BETA vectors of two-electron coupling coefficients for open shells.
  chkpt_wt_iopen(0);

  // Writes out the total energy.
  chkpt_wt_etot(total_energy);
  chkpt_wt_escf(total_energy);
  chkpt_wt_eref(total_energy);


  chkpt_wt_orbspi(sopi);
  chkpt_wt_clsdpi(docc);
  chkpt_wt_openpi(actv);

  int* frz = new int[nirreps];
  for(int h = 0; h < nirreps; ++h) frz[h] = 0;
  chkpt_wt_frzcpi(frz);
  chkpt_wt_frzvpi(frz);
  delete[] frz;

  double** C_save;
  allocate2(double,C_save,nso,nso);
  
  for(int h = 0; h < nirreps; ++h)
    for(int i = 0; i < sopi[h]; ++i)
      for(int j = 0; j < sopi[h]; ++j)
        C_save[i + block_offset[h]][j + block_offset[h]] = C->get(h,i,j);
  chkpt_wt_scf(C_save);

  release2(C_save);

  int k = 0;
  double* evals = new double[nso];
  for(int h = 0; h < nirreps; ++h){
    for(int i = 0; i < sopi[h]; ++i){
      evals[k] = epsilon->get(h,i);
      k++;
    }
  }
  chkpt_wt_evals(evals);
  delete[] evals;
}

}} /* End Namespaces */
