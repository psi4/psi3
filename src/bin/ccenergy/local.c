#include <stdio.h>
#include <math.h>
#include <libciomr.h>
#include <psio.h>
#include <iwl.h>
#include <file30.h>
#include <qt.h>
#include <dpd.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

/* 
** local_init(): Set up parameters of local excitation domains.
**
** Primary parameters include the transformation matrix between the 
** local (projected) virtual space and the orthonormal MO space (V
** matrix), the orbital domain lists, .
*/

void local_init(void)
{
  int i, j, k, ij, stat, a, r, l, I;
  int nmo, nao, nocc, nocc_all, nvir, noei, nirreps, nfzc;
  double **C;  /* AO -> localized MO transformation matrix */
  double **Ci; /* localized MO -> AO transformation matrix */
  double **D;  /* 1/2 SCF closed-shell density matrix (AO) */
  double **Rt; /* Projected, redundant virtual transform (R-tilde) */
  double **S;  /* AO overlap */
  double **St; /* Projected virtual overlap */
  double **Xt; /* Projected, non-redundant virtual transform (X-tilde) */
  double **V;  /* MO -> projected, redundant virtual transform */
  double **Fmo;/* MO basis Fock matrix */
  double **F;  /* AO basis Fock matrix */
  double **Ft; /* Projected, redundant virtual Fock matrix */
  double **Fbar; /* Projected, non-redundant virtual Fock matrix */
  double **W;  /* Transformation matrix from tilde -> bar */
  double *eps; /* orbital energies for local denominators */
  double **X, **Y, **Z;
  double *evals, **evecs;
  double *eps_all; /* All MO energies */
  dpdfile2 fock;

  file30_init();
  C = file30_rd_scf();
  eps_all = file30_rd_evals();
  file30_close();

  /* C1 symmetry only */
  nirreps = moinfo.nirreps;
  if(nirreps != 1) {
    fprintf(outfile, "\nError: localization must use C1 symmetry.\n");
    exit(2);
  }
  /* nao must be nso for now */
  nao = moinfo.nso;
  nmo = moinfo.nmo; /* should be the same as nao */
  if(nmo != nao) {
    fprintf(outfile, "\nError: NMO != NAO!  %d != %d\n", nmo, nao);
    exit(2);
  }
  nocc = moinfo.occpi[0]; /* active doubly occupied orbitals */
  nfzc = moinfo.frdocc[0];  /* frozen doubly occupied orbitals */
  nocc_all = nocc + nfzc; /* all doubly occupied orbitals */
  nvir = moinfo.virtpi[0]; /* active virtual orbitals */

  /* A couple of scratch arrays */
  X = block_matrix(nao, nao);
  Y = block_matrix(nao, nao);

  /* Invert C */
  Ci = block_matrix(nao, nao);
  for(i=0; i < nao; i++)
    for(j=0; j < nao; j++)
      Y[i][j] = C[i][j];

  invert_matrix(C, Ci, nao, outfile);

  for(i=0; i < nao; i++)
    for(j=0; j < nao; j++)
      C[i][j] = Y[i][j];

  fprintf(outfile, "\n\tC inverse (Ci):\n");
  print_mat(Ci, nao, nao, outfile);

  /* Compute the AO-basis overlap integrals */
  S = block_matrix(nao,nao);
  C_DGEMM('t','n',nao,nao,nao,1.0,&(Ci[0][0]),nao,&(Ci[0][0]),nao,
	  0.0,&(S[0][0]),nao);

  fprintf(outfile, "\n\tAO Overlap (S)\n");
  print_mat(S, nao, nao, outfile);

  /* Build the SCF closed-shell density matrix/2 */
  D = block_matrix(nao,nao);
  for(i=0; i < nao; i++) 
    for(j=0; j < nao; j++)
      for(k=0; k < nocc_all; k++)
	D[i][j] += C[i][k] * C[j][k];

  fprintf(outfile, "\n\tAO-basis SCF Density (D):\n");
  print_mat(D, nao, nao, outfile);

  /* Compute the virtual space projector */
  Rt = block_matrix(nao,nao);
  for(i=0; i < nao; i++) Rt[i][i] = 1.0;

  C_DGEMM('n','n',nao,nao,nao,-1.0,&(D[0][0]),nao,&(S[0][0]),nao,
	  1.0,&(Rt[0][0]),nao);

  fprintf(outfile, "\n\tVirtual-Space Projector (R-tilde):\n");
  print_mat(Rt, nao, nao, outfile);

  /* Virtual space metric */
  St = block_matrix(nao,nao);
  C_DGEMM('n','n',nao,nao,nao,1.0,&(S[0][0]),nao,&(Rt[0][0]),nao,
	  0.0,&(X[0][0]),nao);
  C_DGEMM('t','n',nao,nao,nao,1.0,&(Rt[0][0]),nao,&(X[0][0]),nao,
	  0.0,&(St[0][0]),nao);

  fprintf(outfile, "\n\tVirtual-Space Metric (S-tilde):\n");
  print_mat(St, nao, nao, outfile);

  /* Diagonalize metric */
  evals = init_array(nao);
  evecs = block_matrix(nao,nao);
  sq_rsp(nao,nao,St,evals,1,evecs,1e-12);

  fprintf(outfile, "\n\tS-tilde eigenvalues:\n");
  for(i=0; i < nao; i++) fprintf(outfile, "\t%d %20.12f\n", i, evals[i]);

  fprintf(outfile, "\n\tS-tilde eigenvectors:\n");
  print_mat(evecs,nao,nao,outfile);

  /* Build the projected, non-redundant transform (X-tilde) */
  Xt = block_matrix(nao,nao-nocc_all);
  for(i=0,I=0; i < nao; i++) {
    if(evals[i] > 1e-6) {
      for(j=0; j < nao; j++)
	Xt[j][I] = evecs[j][i]/sqrt(evals[i]);
      I++;
    }
  }

  fprintf(outfile, "\n\tTransform to non-redundant, projected virtuals (X-tilde):\n");
  print_mat(Xt, nao, nao-nocc_all, outfile);

  free_block(evecs);
  free(evals);

  V = block_matrix(nvir,nao);
  for(a=0; a < nvir; a++) 
    for(r=0; r < nao; r++) 
      for(i=0; i < nao; i++)
	for(j=0; j < nao; j++)
	  V[a][r] += C[i][a+nocc_all] * S[i][j] * Rt[j][r];

  /* Grab the MO-basis Fock matrix */
  Fmo = block_matrix(nao, nao);
  for(i=0; i < nfzc; i++) Fmo[i][i] = eps_all[i];
  dpd_file2_init(&fock, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_print(&fock, outfile);
  dpd_file2_mat_init(&fock);
  dpd_file2_mat_rd(&fock);
  for(i=0; i < nocc; i++) 
    for(j=0; j < nocc; j++)
      Fmo[i+nfzc][j+nfzc] = fock.matrix[0][i][j];
  dpd_file2_mat_close(&fock);
  dpd_file2_close(&fock);

  dpd_file2_init(&fock, CC_OEI, 0, 1, 1, "fAB");
  dpd_file2_print(&fock, outfile);
  dpd_file2_mat_init(&fock);
  dpd_file2_mat_rd(&fock);
  for(i=0; i < nvir; i++) 
    for(j=0; j < nvir; j++)
      Fmo[i+nfzc+nocc][j+nfzc+nocc] = fock.matrix[0][i][j];
  dpd_file2_mat_close(&fock);
  dpd_file2_close(&fock);

  fprintf(outfile, "\n\tMO Basis Fock matrix:\n");
  print_mat(Fmo, nao, nao, outfile);

  /* Build the AO-basis Fock matrix */
  F = block_matrix(nao,nao);
  C_DGEMM('t','n',nao,nao,nao,1.0,&(Ci[0][0]),nao,&(Fmo[0][0]),nao,
	  0.0,&(X[0][0]),nao);
  C_DGEMM('n','n',nao,nao,nao,1.0,&(X[0][0]),nao,&(Ci[0][0]),nao,
	  0.0,&(F[0][0]),nao);

  fprintf(outfile, "\n\tAO-Basis Fock Matrix:\n");
  print_mat(F, nao, nao, outfile);

  /* Build the projected (redundant) virtual Fock matrix */
  Ft = block_matrix(nao, nao);
  C_DGEMM('t','n',nao,nao,nao,1.0,&(Rt[0][0]),nao,&(F[0][0]),nao,
	  0.0,&(X[0][0]),nao);
  C_DGEMM('n','n',nao,nao,nao,1.0,&(X[0][0]),nao,&(Rt[0][0]),nao,
	  0.0,&(Ft[0][0]),nao);

  /* Project the Fock matrix into the non-redundant virtual space */
  Fbar = block_matrix(nao-nocc_all,nao-nocc_all);
  C_DGEMM('t','n',nao-nocc_all,nao,nao,1.0,&(Xt[0][0]),nao-nocc_all,&(Ft[0][0]),nao,
	  0.0,&(X[0][0]),nao);
  C_DGEMM('n','n',nao-nocc_all,nao-nocc_all,nao,1.0,&(X[0][0]),nao,&(Xt[0][0]),nao-nocc_all,
	  0.0,&(Fbar[0][0]),nao-nocc_all);

  fprintf(outfile, "\n\tFbar matrix:\n");
  print_mat(Fbar,nao-nocc_all,nao-nocc_all,outfile);

  /* Diagonalize Fbar */
  evals = init_array(nao-nocc_all);
  evecs = block_matrix(nao-nocc_all,nao-nocc_all);
  sq_rsp(nao-nocc_all,nao-nocc_all,Fbar,evals,1,evecs,1e-12);

  fprintf(outfile, "\n\tFbar eigenvectors:\n");
  print_mat(evecs,nao-nocc_all,nao-nocc_all,outfile);

  /* Finally, build the W matrix */
  W = block_matrix(nao,nao-nocc_all);
  C_DGEMM('n','n',nao,nao-nocc_all,nao-nocc_all,1.0,&(Xt[0][0]),nao-nocc_all,
	  &(evecs[0][0]),nao-nocc_all,0.0,&(W[0][0]),nao-nocc_all);

  fprintf(outfile, "\n\tW Transformation Matrix:\n");
  print_mat(W,nao,nao-nocc_all,outfile);

  /* build the orbital energy list */
  eps = init_array(nao-nfzc);
  for(i=0; i < nocc; i++)
    eps[i] = Fmo[i+nfzc][i+nfzc];  /* occupied orbital energies */
  for(i=0; i < nao-nocc_all; i++)
    eps[i+nocc] = evals[i]; /* virtual orbital energies */

  fprintf(outfile, "\n\tOrbital Energies:\n");
  for(i=0; i < nao-nfzc; i++)
    fprintf(outfile, "%d %20.12f\n", i, eps[i]);

  free(evals);
  free_block(evecs);


  free_block(F);
  free_block(Fmo);
  free_block(Ft);
  free_block(Fbar);
  free_block(Xt);
  free_block(S);
  free_block(St);
  free_block(Rt);
  free_block(D);
  free_block(Ci);
  free_block(C);

  free_block(X);
  free_block(Y);

  local.V = V;
  local.W = W;
  local.eps = eps;

  fprintf(outfile, "\n\tLocal parameters ready.\n");
}

/*
** local_done(): 
*/

void local_done(void)
{
  free_block(local.V);
  free_block(local.W);
  free(local.eps);

  fprintf(outfile, "\n\tLocal parameters free.\n");
}

void local_filter_T2(dpdbuf4 *T2)
{
  int ij, i, j, a, b;
  int nao, nocc, nvir;
  double **V, **W, *eps;
  double **X, **T2tilde, **T2bar;

  nao = moinfo.nso;
  nocc = moinfo.clsdpi[0];
  nvir = moinfo.virtpi[0];
  V = local.V;
  W = local.W;
  eps = local.eps;

  /* Grab the MO-basis T2's */
  dpd_buf4_mat_irrep_init(T2, 0);
  dpd_buf4_mat_irrep_rd(T2, 0);

  /* Transform the virtuals to the redundant projected virtual basis */
  X = block_matrix(nao,nvir);
  T2tilde = block_matrix(T2->params->rowtot[0], nao*nao);
  for(ij=0; ij < T2->params->rowtot[0]; ij++) {
    C_DGEMM('t', 'n', nao, nvir, nvir, 1.0, &(V[0][0]), nao,
	    &(T2->matrix[0][ij][0]), nvir, 0.0, &(X[0][0]), nvir);
    C_DGEMM('n', 'n', nao, nao, nvir, 1.0, &(X[0][0]), nvir,
	    &(V[0][0]), nao, 0.0, &(T2tilde[ij][0]), nao);
  }
  free_block(X);

  /* Transform the virtuals to the non-redundant virtual basis */
  X = block_matrix(nvir,nao);
  T2bar = block_matrix(T2->params->rowtot[0], nvir*nvir);
  for(ij=0; ij < T2->params->rowtot[0]; ij++) {
    C_DGEMM('t', 'n', nvir, nao, nao, 1.0, &(W[0][0]), nvir,
	    &(T2tilde[ij][0]), nao, 0.0, &(X[0][0]), nao);
    C_DGEMM('n', 'n', nvir, nvir, nao, 1.0, &(X[0][0]), nao,
	    &(W[0][0]), nvir, 0.0, &(T2bar[ij][0]), nvir);
  }
  free_block(X);

  /* Divide the new amplitudes by the denominators */
  for(i=0,ij=0; i < nocc; i++) {
    for(j=0; j < nocc; j++,ij++) {
      for(a=0; a < nvir; a++) {
	for(b=0; b < nvir; b++) {
	  T2bar[ij][a*nvir+b] /= (eps[i] + eps[j] - eps[a+nocc] - eps[b+nocc]);
	}
      }
    }
  }

  /* Transform the new T2's to the redundant virtual basis */
  X = block_matrix(nao, nvir);
  for(ij=0; ij < T2->params->rowtot[0]; ij++) {
    C_DGEMM('n', 'n', nao, nvir, nvir, 1.0, &(W[0][0]), nvir, 
	    &(T2bar[ij][0]), nvir, 0.0, &(X[0][0]), nvir);
    C_DGEMM('n','t', nao, nao, nvir, 1.0, &(X[0][0]), nvir,
	    &(W[0][0]), nvir, 0.0, &(T2tilde[ij][0]), nao);
  }
  free_block(X);
  free_block(T2bar);

  /* Transform the new T2's to the MO basis */
  X = block_matrix(nvir, nao);
  for(ij=0; ij < T2->params->rowtot[0]; ij++) {
    C_DGEMM('n', 'n', nvir, nao, nao, 1.0, &(V[0][0]), nao,
	    &(T2tilde[ij][0]), nao, 0.0, &(X[0][0]), nao);
    C_DGEMM('n', 't', nvir, nvir, nao, 1.0, &(X[0][0]), nao,
	    &(V[0][0]), nao, 0.0, &(T2->matrix[0][ij][0]), nvir);
  }
  free_block(X);
  free_block(T2tilde);

  /* Write the updated MO-basis T2's to disk */
  dpd_buf4_mat_irrep_wrt(T2, 0);
  dpd_buf4_mat_irrep_close(T2, 0);
  /*  dpd_buf4_print(T2, outfile, 1); */
}
