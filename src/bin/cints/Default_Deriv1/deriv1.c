#include<stdio.h>
#include<ip_libv1.h>
#include<math.h>
#include<libciomr.h>
#include<libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include"moinfo.h"
#include"compute_scf_opdm.h"
#include"read_gen_opdm.h"
#include"enuc_deriv1.h"
#include"oe_deriv1.h"
#include"te_deriv1.h"
#include"rot_inv.h"
#include"file11.h"

void deriv1()
{
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
    te_deriv1();
    check_rot_inv();
    if (!strcmp(UserOptions.wfn,"SCF"))
      cleanup_moinfo();
  }

  file11();
  free_block(Grad);

  return;
}

