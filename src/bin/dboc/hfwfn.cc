
#include <stdexcept>
#include <stdlib.h>
extern "C" {
  #include <psifiles.h>
  #include <libchkpt/chkpt.h>
  #include <libciomr/libciomr.h>
}
#include "hfwfn.h"

using namespace std;

HFWavefunction::HFWavefunction()
{
  int unit_opened = 1;
  if (!psio_open_check(PSIF_CHKPT)) {
    chkpt_init(PSIO_OPEN_OLD);
    unit_opened = 0;
  }

  num_mo_ = chkpt_rd_nmo();
  num_so_ = chkpt_rd_nso();
  num_ao_ = chkpt_rd_nao();
  
  refnum_ = (reftype) chkpt_rd_ref();
  
  aotoso_ = chkpt_rd_usotao();
  rref_ = chkpt_rd_rref();
  if (refnum_ == ref_uhf || refnum_ == ref_uks) {
    alpha_evec_ = chkpt_rd_alpha_scf();
    beta_evec_ = chkpt_rd_beta_scf();
  }
  else {
    alpha_evec_ = chkpt_rd_scf();
    beta_evec_ = NULL;
  }

  if (!unit_opened)
    chkpt_close();
}

HFWavefunction::~HFWavefunction()
{
  free_block(aotoso_);
  free_block(rref_);
  free_block(alpha_evec_);
  if (beta_evec_ != NULL) free_block(beta_evec_);
}

int
HFWavefunction::num_ao() { return num_ao_; }

double**
HFWavefunction::alpha_evec() { return alpha_evec_; }

double**
HFWavefunction::beta_evec() {
  if (refnum_ == ref_uhf || refnum_ == ref_uks)
    return beta_evec_;
  else
    return alpha_evec_;
}

double**
HFWavefunction::aotoso() { return aotoso_; }

double**
HFWavefunction::rref() { return rref_; }
