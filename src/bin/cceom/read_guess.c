#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <libipv1/ip_lib.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

struct C1_Guess {
  int i;
  int a;
  double value;
};

void read_guess(int C_irr)
{
  int i, a, k, l, errcod;
  int num_vectors, vector_len, this_irrep, useit;
  char lbl[32];
  double norm, value;
  dpdfile2 CME;
  struct C1_Guess *guess;

  if(params.eom_ref == 2) {
    fprintf(outfile, "\n\tUser-input guesses for UHF states not yet supported.\n");
    exit(2);
  }

  /* User-input guess vectors will override "states_per_irrep" input */
  eom_params.states_per_irrep[C_irr] = 0;

  errcod = ip_count("EOM_GUESS",&num_vectors,0);
  if(errcod != IPE_OK) {
    fprintf(outfile, "\nread_guess(): Unable to read number of guesses from input.\n");
    exit(2);
  }

  for(k=0; k < num_vectors; k++) {

    ip_count("EOM_GUESS", &vector_len, 1, k);
    guess = (struct C1_Guess *) malloc(vector_len * sizeof(struct C1_Guess));

    fprintf(outfile, "\nUser input vector %d:\n", k);
    fprintf(outfile,   "----------------------:\n");
    for(l=0; l < vector_len; l++) {
      errcod = ip_data("EOM_GUESS", "%d", &i, 3, k, l, 0);
      errcod = ip_data("EOM_GUESS", "%d", &a, 3, k, l, 1);
      errcod = ip_data("EOM_GUESS", "%lf", &value, 3, k, l, 2);

      guess[l].i = i;
      guess[l].a = a;
      guess[l].value = value;

      if(l==0) { /* check symmetry of this state */
	this_irrep = moinfo.occ_sym[i]^moinfo.vir_sym[a];
	if(this_irrep != C_irr) {
	  useit = 0;
	  break; /* can't use this guess for this irrep */
	}
	else { 
	  eom_params.states_per_irrep[this_irrep^moinfo.sym]++;
	  useit = 1;
	}
      }

      if(l!=0 && moinfo.occ_sym[i]^moinfo.vir_sym[a] != this_irrep) {
	fprintf(outfile, "\nInconsisent symmetries in components of guess %d.\n", k);
	exit(2);
      }
    }

    if(useit) {
      sprintf(lbl, "%s %d", "CME", k);
      dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
      dpd_file2_scm(&CME, 0);
      dpd_file2_mat_init(&CME);
      norm = 0.0;
      for(l=0; l < vector_len; l++) {
	CME.matrix[C_irr][guess[l].i][guess[l].a] = guess[l].value;
	norm += guess[l].value * guess[l].value;
      }
      norm = sqrt(2.0 * norm);
      for(l=0; l < vector_len; l++)  /* normalize the guess */
	CME.matrix[C_irr][guess[l].i][guess[l].a] *= 1.0/norm;
      dpd_file2_mat_wrt(&CME);
      dpd_file2_mat_close(&CME);
      if(params.eom_ref == 1) {
	sprintf(lbl, "%s %d", "Cme", k);
	dpd_file2_copy(&CME, EOM_Cme, lbl);
      }
      dpd_file2_close(&CME);
    }
    free(guess);
  }

  eom_params.cs_per_irrep[C_irr] = eom_params.states_per_irrep[C_irr];

  fprintf(outfile, "State sought in irrep %d = %d\n", C_irr, eom_params.states_per_irrep[C_irr]);
}
