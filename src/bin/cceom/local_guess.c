#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define EXTERN 
#include "globals.h"

struct onestack {
  int i;
  int a;
  double value;
};

void stack_insert(struct onestack *stack, double value, int i, int a, int level, int stacklen);

void local_guess(void)
{
  int nso, nocc, nvir, nroot, i, ii, a, m;
  int *pairdom_len, *pairdom_nrlen, *weak_pairs;
  double ***V, ***W, *eps_occ, **eps_vir;
  double *T1bar, *T1tilde;
  double fii, value, norm;

  struct onestack *stack;

  char lbl[32];
  dpdfile2 CME;

  nso = local.nso;
  nocc = local.nocc;
  nvir = local.nvir;
  V = local.V;
  W = local.W;
  eps_occ = local.eps_occ;
  eps_vir = local.eps_vir;
  pairdom_len = local.pairdom_len;
  pairdom_nrlen = local.pairdom_nrlen;

  nroot = eom_params.states_per_irrep[0]; /* only C1 allowed */

  stack = (struct onestack *) malloc(nroot * sizeof(struct onestack));
  for(m=0; m < nroot; m++) { stack[m].i = -1; stack[m].a = -1; stack[m].value = 1e12; }

  /* find the nroot lowest excitations in the non-redunant, orthogonal (bar) space */
  for(i=0; i < nocc; i++) {
    ii = i * nocc + i;
    fii = eps_occ[i];
    for(a=0; a < pairdom_nrlen[ii]; a++) {
      value = eps_vir[ii][a] - fii;
      for(m=0; m < nroot; m++) {
	if((fabs(value) < fabs(stack[m].value)) ) {
	  stack_insert(stack, value, i, a, m, nroot);
	  break;
	}
      }
    }
  }

  T1bar = init_array(nso);
  T1tilde = init_array(nso);

  fprintf(outfile, "\n\tTransitions for local guesses:\n");
  fprintf(outfile,   "\t------------------------------\n");
  for(m=0; m < nroot; m++) {
    fprintf(outfile, "\t%3d %3d %14.10f\n", stack[m].i, stack[m].a, stack[m].value);

    memset((void *) T1bar, 0, nso*sizeof(double));
    memset((void *) T1tilde, 0, nso*sizeof(double));

    i = stack[m].i;
    ii = i * nocc + i;

    /* Unit guess vector */
    T1bar[stack[m].a] = 1.0;

    /* Transform this to the canonical MO basis */
    sprintf(lbl, "%s %d", "CME", m);
    dpd_file2_init(&CME, EOM_CME, 0, 0, 1, lbl);
    dpd_file2_mat_init(&CME);

    C_DGEMV('n', pairdom_len[ii], pairdom_nrlen[ii], 1.0, &(W[ii][0][0]), pairdom_nrlen[ii],
	    &(T1bar[0]), 1, 0.0, &(T1tilde[0]), 1);
    C_DGEMV('n', nvir, pairdom_len[ii], 1.0, &(V[ii][0][0]), pairdom_len[ii],
	    &(T1tilde[0]), 1, 0.0, &(CME.matrix[0][i][0]), 1);

    /* normalize this guess in the MO basis */
    norm = 0.0;
    for(a=0; a < nvir; a++) {
      norm += CME.matrix[0][i][a] * CME.matrix[0][i][a];
    }
    norm = sqrt(2.0 * norm);
    fprintf(outfile, "Norm of guess vector %d = %20.14f\n", m, norm);
    for(a=0; a < nvir; a++) {
      CME.matrix[0][i][a] *= 1.0/norm;
    }

    dpd_file2_mat_wrt(&CME);
    dpd_file2_mat_close(&CME);

    dpd_file2_close(&CME);
  }

  fprintf(outfile, "\n");

  free(T1bar);
  free(T1tilde);

  free(stack);

  eom_params.cs_per_irrep[0] = nroot;
}

void stack_insert(struct onestack *stack, double value, int i, int a, int level, int stacklen)
{
  int l;
  struct onestack temp;

  temp = stack[level];

  stack[level].value = value;
  stack[level].i = i;
  stack[level].a = a;

  value = temp.value;
  i = temp.i;
  a = temp.a;

  for(l=level; l < stacklen-1; l++) {
    temp = stack[l+1];

    stack[l+1].value = value;
    stack[l+1].i = i;
    stack[l+1].a = a;

    value = temp.value;
    i = temp.i;
    a = temp.a;
  }
}

