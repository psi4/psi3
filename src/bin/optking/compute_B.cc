// This function constructs the B matrix for a set of salcs

extern "C" {
  #include <stdio.h>
  #include <file30.h>
  #include <stdlib.h>
  #include <string.h>
  #include <physconst.h>
  #include <math.h>
  #include <libciomr.h>
}

#define EXTERN
#include "opt.h"
#include "cartesians.h"
#include "internals.h"
#include "salc.h"

double **compute_B(int num_atoms,internals &simples,salc_set &symm) {
  int i,j,k,a,b,c,d, simple, intco_type, sub_index;
  double **B, coeff, prefactor;

  B = init_matrix(symm.get_num(),num_atoms*3);

  for (i=0;i<symm.get_num();++i) {
     prefactor = symm.get_prefactor(i);
     for (j=0;j<symm.get_length(i);++j) {
        simple = symm.get_simple(i,j);
        coeff = symm.get_coeff(i,j);
        simples.locate_id(simple,&intco_type,&sub_index);
        if (intco_type == STRE_TYPE) {
           a = simples.stre.get_A(sub_index);
           b = simples.stre.get_B(sub_index);
           for (k=0;k<3;++k) {
              B[i][3*a+k] += prefactor * coeff * simples.stre.get_s_A(sub_index,k);
              B[i][3*b+k] += prefactor * coeff * simples.stre.get_s_B(sub_index,k);
           }
        }
        else if (intco_type == BEND_TYPE) {
           a = simples.bend.get_A(sub_index);
           b = simples.bend.get_B(sub_index);
           c = simples.bend.get_C(sub_index);
           for (k=0;k<3;++k) {
              B[i][3*a+k] += prefactor * coeff * simples.bend.get_s_A(sub_index,k);
              B[i][3*b+k] += prefactor * coeff * simples.bend.get_s_B(sub_index,k);
              B[i][3*c+k] += prefactor * coeff * simples.bend.get_s_C(sub_index,k);
           }
        }
        else if (intco_type == TORS_TYPE) {
           a = simples.tors.get_A(sub_index);
           b = simples.tors.get_B(sub_index);
           c = simples.tors.get_C(sub_index);
           d = simples.tors.get_D(sub_index);
           for (k=0;k<3;++k) {
              B[i][3*a+k] += prefactor * coeff * simples.tors.get_s_A(sub_index,k);
              B[i][3*b+k] += prefactor * coeff * simples.tors.get_s_B(sub_index,k);
              B[i][3*c+k] += prefactor * coeff * simples.tors.get_s_C(sub_index,k);
              B[i][3*d+k] += prefactor * coeff * simples.tors.get_s_D(sub_index,k);
           }
        }
        else if (intco_type == OUT_TYPE) {
           a = simples.out.get_A(sub_index);
           b = simples.out.get_B(sub_index);
           c = simples.out.get_C(sub_index);
           d = simples.out.get_D(sub_index);
           for (k=0;k<3;++k) {
              B[i][3*a+k] += prefactor * coeff * simples.out.get_s_A(sub_index,k);
              B[i][3*b+k] += prefactor * coeff * simples.out.get_s_B(sub_index,k);
              B[i][3*c+k] += prefactor * coeff * simples.out.get_s_C(sub_index,k);
              B[i][3*d+k] += prefactor * coeff * simples.out.get_s_D(sub_index,k);
           }
        }
     }
  }

  return B;
}

