/*! \file 
    \ingroup (CINTS)
    \brief Enter brief description of file here 
*/
#include<stdio.h>
#include<libipv1/ip_lib.h>
#include<math.h>
#include<libciomr/libciomr.h>
#include<libint/libint.h>
#include<pthread.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include"moinfo.h"
#include"compute_scf_opdm.h"
#include"read_gen_opdm.h"
#include"enuc_deriv1.h"
#include"oe_deriv1.h"
#include"te_deriv1_scf.h"
#include"te_deriv1_corr.h"
#include"symmetrize_deriv1.h"
#include"rot_inv.h"
#include"file11.h"

pthread_mutex_t deriv1_mutex;
double **grad_te;

void deriv1()
{
  /* Either contract integrals with the densities ... */
  if (UserOptions.symm_ints == 0) {
    /*--- Gradient in the canonical frame ---*/
    Grad = block_matrix(Molecule.num_atoms,3);
    
    if (Molecule.num_atoms != 0) {
      if (!strcmp(UserOptions.wfn,"SCF")) {
        init_moinfo();
        compute_scf_opdm();
      }
      else
        read_gen_opdm();
      enuc_deriv1();
      oe_deriv1();
      if (!strcmp(UserOptions.wfn,"SCF"))
        te_deriv1_scf();
      else
        te_deriv1_corr();
      symmetrize_deriv1();
      check_rot_inv();
      if (!strcmp(UserOptions.wfn,"SCF"))
        cleanup_moinfo();
    }
    
    file11();
    chkpt_wt_grad(Grad[0]);
    free_block(Grad);
  }
  /* ... or compute integrals alone */
  else {
    oe_deriv1_ints();
    te_deriv1_ints();
  }

  return;
}

