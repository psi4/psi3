#include <stdio.h>
#include <math.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libdpd/dpd.h>
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

double **fock_build(double **);

/* rotate(): Rotate the orbitals using a linear transformation matrix
** built from converged T1 amplitudes.  I still need to add spin-restricted
** Brueckner rotations from my 1997 paper.
**
** TDC, 5/03
*/

int rotate(void)
{
  int i, a, ii, aa, j, ij;
  int p, q;
  int h, nirreps, nso, nmo, ntri, stat;
  dpdfile2 T1;
  double **U, **S, **X, *scratch;
  double *evals, *work, **SO_S, **MO_S;
  double **scf, **scf_new, **scf_final;
  double max;
  double **D; /* SCF density */
  double **fock; /* Fock matrix (SO or MO basis) */
  int *offset;

  nirreps = moinfo.nirreps;
  nso = moinfo.nso;
  nmo = moinfo.nmo;

  /* First check to see if we've already converged the orbitals */
  max = 0.0;
  if(params.ref == 0) { /** RHF **/
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_mat_init(&T1);
    dpd_file2_mat_rd(&T1);

    for(h=0; h < nirreps; h++)
      for(i=0; i < moinfo.occpi[h]; i++)
	for(a=0; a < moinfo.virtpi[h]; a++)
	  if(fabs(T1.matrix[h][i][a]) > max) max = T1.matrix[h][i][a];

    dpd_file2_mat_close(&T1);
    dpd_file2_close(&T1);
  }
  else if(params.ref == 1) { /** UHF **/

    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_mat_init(&T1);
    dpd_file2_mat_rd(&T1);

    for(h=0; h < nirreps; h++)
      for(i=0; i < moinfo.aoccpi[h]; i++)
	for(a=0; a < moinfo.avirtpi[h]; a++)
	  if(fabs(T1.matrix[h][i][a]) > max) max = T1.matrix[h][i][a];

    dpd_file2_mat_close(&T1);
    dpd_file2_close(&T1);

    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_file2_mat_init(&T1);
    dpd_file2_mat_rd(&T1);

    for(h=0; h < nirreps; h++)
      for(i=0; i < moinfo.boccpi[h]; i++)
	for(a=0; a < moinfo.bvirtpi[h]; a++)
	  if(fabs(T1.matrix[h][i][a]) > max) max = T1.matrix[h][i][a];

    dpd_file2_mat_close(&T1);
    dpd_file2_close(&T1);
  }

  if(fabs(max) <= params.bconv) {
    fprintf(outfile, "\tBrueckner orbitals converged.\n");
    return(1);
  }
  else 
    fprintf(outfile, "\tRotating orbitals.  Maximum T1 = %15.12f\n", max);

  /* grab the SO-basis overlap integrals for later use */
  SO_S = block_matrix(nso, nso);
  ntri = nso * (nso+1)/2;
  scratch = init_array(ntri);
  stat = iwl_rdone(PSIF_OEI, PSIF_SO_S, scratch, ntri, 0, 0, outfile);
  for(i=0,ij=0; i < nso; i++)
    for(j=0; j <= i; j++,ij++) {
      SO_S[i][j] = SO_S[j][i] = scratch[ij];
    }
  free(scratch);

  if(params.ref == 0) { /* RHF */

    U = block_matrix(nmo, nmo);
    for(i=0; i < nmo; i++) U[i][i] = 1.0;

    max = 0.0;
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_mat_init(&T1);
    dpd_file2_mat_rd(&T1);
    for(h=0; h < nirreps; h++) {
      for(i=0; i < moinfo.occpi[h]; i++) {
	ii = moinfo.qt2pitzer[moinfo.qt_occ[i] + moinfo.occ_off[h]];
	for(a=0; a < moinfo.virtpi[h]; a++) {
	  aa = moinfo.qt2pitzer[moinfo.qt_vir[a] + moinfo.vir_off[h]];

	  U[ii][aa] = T1.matrix[h][i][a];
	  U[aa][ii] = -T1.matrix[h][i][a];
	}
      }
    }
    dpd_file2_mat_close(&T1);
    dpd_file2_close(&T1);

    chkpt_init(PSIO_OPEN_OLD);
    scf = chkpt_rd_scf();

    /* build the SO-basis SCF density */
    offset = init_int_array(nirreps);
    for(h=1; h < nirreps; h++)
      offset[h] = offset[h-1] + moinfo.orbspi[h-1];

    D = block_matrix(nso,nso);
    for(h=0; h < nirreps; h++)
      for(p=offset[h]; p < offset[h]+moinfo.orbspi[h]; p++)
	for(q=offset[h]; q < offset[h]+moinfo.orbspi[h]; q++)
	  for(i=offset[h]; i < offset[h]+moinfo.occpi[h]; i++)
	    D[p][q] += scf[p][i] * scf[q][i];
    free(offset);

    /* build the SO-basis Fock matrix */
    fock = fock_build(D);
    free_block(D);

    /* transform the fock matrix to the MO basis */
    X = block_matrix(nso,nso);
    C_DGEMM('n','n',nso,nmo,nso,1.0,&(fock[0][0]),nso,&(scf[0][0]),nmo,
	0,&(X[0][0]),nso);
    C_DGEMM('t','n',nmo,nmo,nso,1.0,&(scf[0][0]),nmo,&(X[0][0]),nso,
	0,&(fock[0][0]),nso);
    free_block(X);
    
    scf_new = block_matrix(nso, nmo);
    C_DGEMM('n','t',nso,nmo,nmo,1,&(scf[0][0]),nmo,&(U[0][0]),nmo,
	    0,&(scf_new[0][0]),nmo);

    MO_S = block_matrix(nmo, nmo);

    X = block_matrix(nso, nso);
    C_DGEMM('t','n',nmo, nso, nso, 1, &(scf_new[0][0]), nmo, &(SO_S[0][0]), nso,
	    0, &(X[0][0]), nso);
    C_DGEMM('n','n',nmo, nmo, nso, 1, &(X[0][0]), nso, &(scf_new[0][0]), nmo, 
	    0, &(MO_S[0][0]), nmo);
    free_block(X);

    evals = init_array(nmo);
    work = init_array(nmo*3);
    if(stat = C_DSYEV('v','u', nmo,&(MO_S[0][0]),nmo,evals,work,nmo*3)) {
      fprintf(outfile, "rotate(): Error in overlap diagonalization. stat = %d\n", stat);
      exit(PSI_RETURN_FAILURE);
    }

    /* build S^-1/2 for this basis */
    S = block_matrix(nmo, nmo);
    for(i=0; i < nmo; i++) {
      if(fabs(evals[i]) > 1e-8) S[i][i] = 1/sqrt(evals[i]);
      else S[i][i] = 0.0;
    }

    X = block_matrix(nmo, nmo);
    C_DGEMM('t','n',nmo, nmo, nmo, 1, &(MO_S[0][0]), nso, &(S[0][0]), nmo,
	    0, &(X[0][0]), nmo);
    C_DGEMM('n','n', nmo, nmo, nmo, 1, &(X[0][0]), nmo, &(MO_S[0][0]), nso, 
	    0, &(S[0][0]), nmo);
    free_block(X);

    /* orthogonalize the basis */
    scf_final = block_matrix(nso, nmo);
    C_DGEMM('n','n',nmo,nmo,nmo,1,&(scf_new[0][0]),nmo,&(S[0][0]),nmo,
	    0,&(scf_final[0][0]),nmo);
    free_block(S);

    free_block(MO_S);

    chkpt_wt_scf(scf_final);
    chkpt_close();

    free_block(scf);
    free_block(scf_new);
    free_block(scf_final);

  }
  else if(params.ref == 2) { /* UHF */

    /* AA block */
    U = block_matrix(nmo, nmo);

    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_mat_init(&T1);
    dpd_file2_mat_rd(&T1);
    for(h=0; h < nirreps; h++) {
      for(i=0; i < moinfo.aoccpi[h]; i++) {
	ii = moinfo.qt2pitzer[moinfo.qt_aocc[i] + moinfo.aocc_off[h]];
	for(a=0; a < moinfo.avirtpi[h]; a++) {
	  aa = moinfo.qt2pitzer[moinfo.qt_avir[a] + moinfo.avir_off[h]];

	  U[ii][aa] = T1.matrix[h][i][a];
	  U[aa][ii] = -T1.matrix[h][i][a];
	}
      }
    }
    dpd_file2_mat_close(&T1);
    dpd_file2_close(&T1);

    chkpt_init(PSIO_OPEN_OLD);
    scf = chkpt_rd_alpha_scf();

    scf_new = block_matrix(nso, nmo);
    C_DGEMM('n','t',nso,nmo,nmo,1,&(scf[0][0]),nmo,&(U[0][0]),nmo,
	    0,&(scf_new[0][0]),nmo);

    MO_S = block_matrix(nmo, nmo);

    X = block_matrix(nso, nso);
    C_DGEMM('t','n',nmo, nso, nso, 1, &(scf_new[0][0]), nmo, &(SO_S[0][0]), nso,
	    0, &(X[0][0]), nso);
    C_DGEMM('n','n',nmo, nmo, nso, 1, &(X[0][0]), nso, &(scf_new[0][0]), nmo, 
	    0, &(MO_S[0][0]), nmo);
    free_block(X);

    evals = init_array(nmo);
    work = init_array(nmo*3);
    if(stat = C_DSYEV('v','u', nmo,&(MO_S[0][0]),nmo,evals,work,nmo*3)) {
      fprintf(outfile, "rotate(): Error in overlap diagonalization. stat = %d\n", stat);
      exit(PSI_RETURN_FAILURE);
    }

    /* build S^-1/2 for this basis */
    S = block_matrix(nmo, nmo);
    for(i=0; i < nmo; i++) {
      if(fabs(evals[i]) > 1e-8) S[i][i] = 1/sqrt(evals[i]);
      else S[i][i] = 0.0;
    }

    X = block_matrix(nmo, nmo);
    C_DGEMM('t','n',nmo, nmo, nmo, 1, &(MO_S[0][0]), nso, &(S[0][0]), nmo,
	    0, &(X[0][0]), nmo);
    C_DGEMM('n','n', nmo, nmo, nmo, 1, &(X[0][0]), nmo, &(MO_S[0][0]), nso, 
	    0, &(S[0][0]), nmo);
    free_block(X);

    /* orthogonalize the basis */
    scf_final = block_matrix(nso, nmo);
    C_DGEMM('n','n',nmo,nmo,nmo,1,&(scf_new[0][0]),nmo,&(S[0][0]),nmo,
	    0,&(scf_final[0][0]),nmo);
    free_block(S);

    free_block(MO_S);

    chkpt_wt_beta_scf(scf_final);
    chkpt_close();

    free_block(scf);
    free_block(scf_new);
    free_block(scf_final);

    /* BB block */
    U = block_matrix(nmo, nmo);

    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_file2_mat_init(&T1);
    dpd_file2_mat_rd(&T1);
    for(h=0; h < nirreps; h++) {
      for(i=0; i < moinfo.boccpi[h]; i++) {
	ii = moinfo.qt2pitzer[moinfo.qt_bocc[i] + moinfo.bocc_off[h]];
	for(a=0; a < moinfo.bvirtpi[h]; a++) {
	  aa = moinfo.qt2pitzer[moinfo.qt_bvir[a] + moinfo.bvir_off[h]];

	  U[ii][aa] = T1.matrix[h][i][a];
	  U[aa][ii] = -T1.matrix[h][i][a];
	}
      }
    }
    dpd_file2_mat_close(&T1);
    dpd_file2_close(&T1);

    chkpt_init(PSIO_OPEN_OLD);
    scf = chkpt_rd_beta_scf();

    scf_new = block_matrix(nso, nmo);
    C_DGEMM('n','t',nso,nmo,nmo,1,&(scf[0][0]),nmo,&(U[0][0]),nmo,
	    0,&(scf_new[0][0]),nmo);

    MO_S = block_matrix(nmo, nmo);

    X = block_matrix(nso, nso);
    C_DGEMM('t','n',nmo, nso, nso, 1, &(scf_new[0][0]), nmo, &(SO_S[0][0]), nso,
	    0, &(X[0][0]), nso);
    C_DGEMM('n','n',nmo, nmo, nso, 1, &(X[0][0]), nso, &(scf_new[0][0]), nmo, 
	    0, &(MO_S[0][0]), nmo);
    free_block(X);

    evals = init_array(nmo);
    work = init_array(nmo*3);
    if(stat = C_DSYEV('v','u', nmo,&(MO_S[0][0]),nmo,evals,work,nmo*3)) {
      fprintf(outfile, "rotate(): Error in overlap diagonalization. stat = %d\n", stat);
      exit(PSI_RETURN_FAILURE);
    }

    /* build S^-1/2 for this basis */
    S = block_matrix(nmo, nmo);
    for(i=0; i < nmo; i++) {
      if(fabs(evals[i]) > 1e-8) S[i][i] = 1/sqrt(evals[i]);
      else S[i][i] = 0.0;
    }

    X = block_matrix(nmo, nmo);
    C_DGEMM('t','n',nmo, nmo, nmo, 1, &(MO_S[0][0]), nso, &(S[0][0]), nmo,
	    0, &(X[0][0]), nmo);
    C_DGEMM('n','n', nmo, nmo, nmo, 1, &(X[0][0]), nmo, &(MO_S[0][0]), nso, 
	    0, &(S[0][0]), nmo);
    free_block(X);

    /* orthogonalize the basis */
    scf_final = block_matrix(nso, nmo);
    C_DGEMM('n','n',nmo,nmo,nmo,1,&(scf_new[0][0]),nmo,&(S[0][0]),nmo,
	    0,&(scf_final[0][0]),nmo);
    free_block(S);

    free_block(MO_S);

    chkpt_wt_beta_scf(scf_final);
    chkpt_close();

    free_block(scf);
    free_block(scf_new);
    free_block(scf_final);
  }

  return 0;
}
