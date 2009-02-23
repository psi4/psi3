#include <libchkpt/chkpt.h>
#include <libmoinfo/libmoinfo.h>

#include "scf.h"

namespace psi{ namespace mcscf{

void SCF::save_info()
{
  // Writes out the total number of irreps in the point group in which the molecule is being considered which have non-zero number of basis functions.
  int n_so_typs = 0;
  for(int h = 0; h < nirreps; ++h){
    if( docc[h] + actv[h] > 0 ) n_so_typs++;
  }
  chkpt_wt_nsymhf(n_so_typs);

  // Write out the total number of molecular orbitals.
  chkpt_wt_nmo(nso); // TODO: find nmo

  // Write out the dimensionality of ALPHA and BETA vectors of two-electron coupling coefficients for open shells.
  int tmp_iopen = ioff[moinfo_scf->get_nactv()];
  if(reference == tcscf) tmp_iopen = -tmp_iopen;
  chkpt_wt_iopen(tmp_iopen);

  // Write open-shell coupling coefficients
  if(moinfo_scf->get_nactv() > 0){
    double** ccvecs;
    allocate2(double,ccvecs,2,ioff[moinfo_scf->get_nactv()]);
    for(int i=0; i < ioff[reference == tcscf]; i++) {
      ccvecs[0][i] = 0.0;
      ccvecs[1][i] = 0.0;
    }
    chkpt_wt_ccvecs(ccvecs);
    release2(ccvecs);
  }

  // Writes out the total energy.
  chkpt_wt_etot(total_energy);
  chkpt_wt_escf(total_energy);
  chkpt_wt_eref(total_energy);


  chkpt_wt_orbspi(&sopi[0]);
  chkpt_wt_clsdpi(&docc[0]);
  chkpt_wt_openpi(&actv[0]);



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
