
// This function will return the irrep of a SALC

extern "C" {
  #include <stdio.h>
  #include <file30.h>
  #include <stdlib.h>
  #include <string.h>
  #include <physconst.h>
  #include <math.h>
  #include <libciomr.h>
  #include <ip_libv1.h>
}

#define EXTERN
#include "opt.h"
#include "cartesians.h"
#include "internals.h"
#include "salc.h"

int irrep(internals &simples, double *salc) {
  int i, j, a, ops, **ct_tors, **ict, num_irreps;
  double sum, *tmp_evect, *evect_character;
  char buffer[MAX_LINELENGTH], *err;
  int h, character_is_set;

  num_irreps = syminfo.num_irreps;
//  evect_character = init_array(num_irreps);
evect_character = (double *) malloc(num_irreps*sizeof(double));
  tmp_evect = init_array(simples.get_num());

  /* Determine the characters of the salc under different operations */
  for (ops = 0;ops<num_irreps;++ops) {
    zero_arr(tmp_evect, simples.get_num());
    for (i=0;i<simples.get_num();++i) {
      a = simples.id_to_index(syminfo.ict_ops[i][ops]);
      tmp_evect[a] += syminfo.ict_ops_sign[i][ops] * salc[i];
    }
//fprintf(outfile,"tmp_evect: for operation %d\n",ops);
//for (i=0;i<simples.get_num();++i) fprintf(outfile,"%15.10lf\n",tmp_evect[i]);
    character_is_set = 0;
    for (i=0;i<simples.get_num();++i) {
      if (fabs(salc[i]) > optinfo.irrep_tol) {
        if (!character_is_set) {
           evect_character[ops] = tmp_evect[i]/salc[i];
           character_is_set = 1;
        }
        else if ((fabs(evect_character[ops]-tmp_evect[i]/salc[i])>optinfo.irrep_tol)) {
          fprintf(outfile,"After character_is_set to %15.10lf,\n",
                  evect_character[ops]);
          fprintf(outfile,"tmp_evect[%d]/salc[%d]: %15.10lf\n", 
                  i, i, tmp_evect[i]/salc[i]);
          exit(2);
        }
      }
    }
  }

  /* Determine the irrep of the salc */
  for (i=0;i<num_irreps;++i)
    for (j=0;j<num_irreps;++j) {
      if (fabs(evect_character[j] - syminfo.ct[i][j]) > optinfo.irrep_tol) break;
      if (j == (num_irreps-1)) {
        free(evect_character);
        free(tmp_evect);
//fprintf(outfile,"returnings irrep: %d\n",i);
//fflush(outfile);
        return i;
      }
    }

  fprintf(outfile,"Eigenvector is not an irrep of point group.\n");
  exit(1);
  return 9999;
}


