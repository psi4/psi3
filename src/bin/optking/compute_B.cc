/*! \file
    \ingroup OPTKING
    \brief This function constructs the B matrix for a set of salcs
*/

#define EXTERN
#include "globals.h"
#undef EXTERN
#include "cartesians.h"
#include "simples.h"
#include "salc.h"
#include "opt.h"

namespace psi { namespace optking {

double **compute_B(const simples_class & simples, const salc_set & symm) {
  int i,j,k,a,b,c,d, simple, sub_index, sub_index2;
  int J,K,atom,xyz,cnt_frag=-1;
  double **B, coeff, prefactor, weight;
  Intco_type intco_type;

  B = block_matrix(symm.get_num(),3*optinfo.natom);

  for (i=0;i<symm.get_num();++i) {
    prefactor = symm.get_prefactor(i);
    for (j=0;j<symm.get_length(i);++j) {
      simple = symm.get_simple(i,j);
      coeff = symm.get_coeff(i,j);
      simples.locate_id(simple,&intco_type,&sub_index,&sub_index2);
      if (intco_type == STRE_TYPE) {
        a = simples.stre[sub_index].get_A();
        b = simples.stre[sub_index].get_B();
        for (k=0;k<3;++k) {
          B[i][3*a+k] += prefactor * coeff * simples.stre[sub_index].get_s_A(k);
          B[i][3*b+k] += prefactor * coeff * simples.stre[sub_index].get_s_B(k);
        }
      }
      else if (intco_type == BEND_TYPE) {
        a = simples.bend[sub_index].get_A();
        b = simples.bend[sub_index].get_B();
        c = simples.bend[sub_index].get_C();
        for (k=0;k<3;++k) {
          B[i][3*a+k] += prefactor * coeff * simples.bend[sub_index].get_s_A(k);
          B[i][3*b+k] += prefactor * coeff * simples.bend[sub_index].get_s_B(k);
          B[i][3*c+k] += prefactor * coeff * simples.bend[sub_index].get_s_C(k);
        }
      }
      else if (intco_type == TORS_TYPE) {
        a = simples.tors[sub_index].get_A();
        b = simples.tors[sub_index].get_B();
        c = simples.tors[sub_index].get_C();
        d = simples.tors[sub_index].get_D();
        for (k=0;k<3;++k) {
          B[i][3*a+k] += prefactor * coeff * simples.tors[sub_index].get_s_A(k);
          B[i][3*b+k] += prefactor * coeff * simples.tors[sub_index].get_s_B(k);
          B[i][3*c+k] += prefactor * coeff * simples.tors[sub_index].get_s_C(k);
          B[i][3*d+k] += prefactor * coeff * simples.tors[sub_index].get_s_D(k);
        }
      }
      else if (intco_type == OUT_TYPE) {
        a = simples.out[sub_index].get_A();
        b = simples.out[sub_index].get_B();
        c = simples.out[sub_index].get_C();
        d = simples.out[sub_index].get_D();
        for (k=0;k<3;++k) {
          B[i][3*a+k] += prefactor * coeff * simples.out[sub_index].get_s_A(k);
          B[i][3*b+k] += prefactor * coeff * simples.out[sub_index].get_s_B(k);
          B[i][3*c+k] += prefactor * coeff * simples.out[sub_index].get_s_C(k);
          B[i][3*d+k] += prefactor * coeff * simples.out[sub_index].get_s_D(k);
        }
      }
      else if (intco_type == LINB_TYPE) {
        a = simples.linb[sub_index].get_A();
        b = simples.linb[sub_index].get_B();
        c = simples.linb[sub_index].get_C();
        for (k=0;k<3;++k) {
          B[i][3*a+k] += prefactor * coeff * simples.linb[sub_index].get_s_A(k);
          B[i][3*b+k] += prefactor * coeff * simples.linb[sub_index].get_s_B(k);
          B[i][3*c+k] += prefactor * coeff * simples.linb[sub_index].get_s_C(k);
        }
      }
      else if (intco_type == FRAG_TYPE) {
        for (a=0; a<simples.frag[sub_index].get_A_natom(); ++a) { /* loop over atoms in A */
          atom   = simples.frag[sub_index].get_A_atom(a);        /* atom number of a'th atom in A */
          for (K=0;K<simples.frag[sub_index].get_A_P();++K) {     /* loop over reference atoms of A */
            weight = simples.frag[sub_index].get_A_weight(K,a);  /* weight of a'th atom in A for this K */
            for (xyz=0;xyz<3;++xyz)
              B[i][3*atom+xyz] += prefactor * coeff * weight *
                simples.frag[sub_index].get_A_s(sub_index2,3*K+xyz);
          }
        }
        for (b=0; b<simples.frag[sub_index].get_B_natom(); ++b) { /* loop over atoms in B */
          atom   = simples.frag[sub_index].get_B_atom(b);        /* atom number of b'th atom in B */
          for (K=0;K<simples.frag[sub_index].get_B_P();++K) { /* loop over reference atoms of B*/
            weight = simples.frag[sub_index].get_B_weight(K,b);  /* weight of b'th atom in B for this K */
            for (xyz=0;xyz<3;++xyz)
              B[i][3*atom+xyz] += prefactor * coeff * weight *
                simples.frag[sub_index].get_B_s(sub_index2,3*K+xyz);
          }
        }
      }
    }
  }

//  fprintf(outfile, "B matrix\n");
//  print_mat2(B, symm.get_num(),3*optinfo.natom, outfile);

  return B;
}

}} /* namespace psi::optking */
