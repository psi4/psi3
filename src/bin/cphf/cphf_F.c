#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

/* cphf_F(): Solve the first-order CPHF equations for an electric
** field perturbation.
**
** The CPHF equations are:
**
** A U^a = B0^a
**
** where A is the MO hessian [computed in mohess()], U^a is the
** orbital response (CPHF) coefficient, and B0^a is the
** perturbation-dependent inhomogenous factor.  For electric field
** perturbations, B0^a is given by (MO basis):
**
** (B0^a)_ai = (mu^a)_ai
**
** where (mu^a)_pq is a dipole moment integral.
**
** For more details (and clearer notation) see:
**
** Y. Yamaguchi et al., "A New Dimension to Quantum Chemistry:
** Analytic Derivative Methods in Ab Initio Molecular Electronic
** Structure Theory", Oxford Press, New York, 1994.  Ch.10,
** pp. 128-132 and Ch.16, pp. 332-333.
**
** TDC, October 2002
*/

void cphf_F(double **A, double ***U)
{
  int noei, stat, coord;
  double ***B, *scratch, **B0, **Acopy;
  double **TMP, **X;
  char *label;
  int a, asym, afirst, alast;
  int i, isym, ifirst, ilast;
  int AI, j, ij;
  int *ipiv, error;

  /* Allocate space for the B vectors */
  B = (double ***) malloc(3 * sizeof(double **));
  for(coord=0; coord < 3; coord++) B[coord] = block_matrix(nmo,nmo);

  /* Transform dipole integrals to the MO basis and save on disk */
  TMP = block_matrix(nao, nao);
  X = block_matrix(nao, nao);
  scratch = init_array(noei_ao);

  stat = iwl_rdone(PSIF_OEI, PSIF_AO_MX, scratch, noei_ao, 0, 0, outfile);
  for(i=0,ij=0; i < nao; i++)
    for(j=0; j <= i; j++,ij++) {
      TMP[i][j] = TMP[j][i] = scratch[ij];
    }

  /*
  fprintf(outfile, "\tAO-basis MuX Integrals:\n");
  print_mat(TMP, nao, nao, outfile);
  */

  C_DGEMM('n','t',nao,nso,nao,1,&(TMP[0][0]),nao,&(usotao[0][0]),nao,
	  0,&(X[0][0]),nao);
  C_DGEMM('n','n',nso,nso,nao,1,&(usotao[0][0]),nao,&(X[0][0]),nao,
	  0,&(TMP[0][0]),nao);

  C_DGEMM('n','n',nso,nmo,nso,1,&(TMP[0][0]),nao,&(scf[0][0]),nmo,
	  0,&(X[0][0]),nao);
  C_DGEMM('t','n',nmo,nmo,nso,1,&(scf[0][0]),nmo,&(X[0][0]),nao,
	  0,&(B[0][0][0]),nmo);

  for(i=0,ij=0; i < nmo; i++)
    for(j=0; j <= i; j++,ij++)
      scratch[ij] = B[0][i][j];
  iwl_wrtone(PSIF_OEI, PSIF_MO_MX, ntri, scratch);

  zero_arr(scratch,noei_ao);
  stat = iwl_rdone(PSIF_OEI, PSIF_AO_MY, scratch, noei_ao, 0, 0, outfile);
  for(i=0,ij=0; i < nao; i++)
    for(j=0; j <= i; j++,ij++) {
      TMP[i][j] = TMP[j][i] = scratch[ij];
    }

  /*
  fprintf(outfile, "\tAO-basis MuY Integrals:\n");
  print_mat(TMP, nao, nao, outfile);
  */

  C_DGEMM('n','t',nao,nso,nao,1,&(TMP[0][0]),nao,&(usotao[0][0]),nao,
	  0,&(X[0][0]),nao);
  C_DGEMM('n','n',nso,nso,nao,1,&(usotao[0][0]),nao,&(X[0][0]),nso,
	  0,&(TMP[0][0]),nao);

  C_DGEMM('n','n',nso,nmo,nso,1,&(TMP[0][0]),nao,&(scf[0][0]),nmo,
	  0,&(X[0][0]),nao);
  C_DGEMM('t','n',nmo,nmo,nso,1,&(scf[0][0]),nmo,&(X[0][0]),nao,
	  0,&(B[1][0][0]),nmo);

  for(i=0,ij=0; i < nmo; i++)
    for(j=0; j <= i; j++,ij++)
      scratch[ij] = B[1][i][j];
  iwl_wrtone(PSIF_OEI, PSIF_MO_MY, ntri, scratch);

  zero_arr(scratch,noei_ao);
  stat = iwl_rdone(PSIF_OEI, PSIF_AO_MZ, scratch, noei_ao, 0, 0, outfile);
  for(i=0,ij=0; i < nao; i++)
    for(j=0; j <= i; j++,ij++) {
      TMP[i][j] = TMP[j][i] = scratch[ij];
    }

  /*
  fprintf(outfile, "\tAO-basis MuZ Integrals:\n");
  print_mat(TMP, nao, nao, outfile);
  */

  C_DGEMM('n','t',nao,nso,nao,1,&(TMP[0][0]),nao,&(usotao[0][0]),nao,
	  0,&(X[0][0]),nao);
  C_DGEMM('n','n',nso,nso,nao,1,&(usotao[0][0]),nao,&(X[0][0]),nso,
	  0,&(TMP[0][0]),nao);

  C_DGEMM('n','n',nso,nmo,nso,1,&(TMP[0][0]),nao,&(scf[0][0]),nmo,
	  0,&(X[0][0]),nao);
  C_DGEMM('t','n',nmo,nmo,nso,1,&(scf[0][0]),nmo,&(X[0][0]),nao,
	  0,&(B[2][0][0]),nmo);

  for(i=0,ij=0; i < nmo; i++)
    for(j=0; j <= i; j++,ij++)
      scratch[ij] = B[2][i][j];
  iwl_wrtone(PSIF_OEI, PSIF_MO_MZ, ntri, scratch);

  free_block(TMP);
  free_block(X);
  free(scratch);

  /*
  fprintf(outfile, "\tSCF MO's:\n");
  print_mat(scf, nso, nmo, outfile);
  fprintf(outfile, "\tMO-basis MuX Integrals:\n");
  print_mat(B[0], nmo, nmo, outfile);
  fprintf(outfile, "\tMO-basis MuY Integrals:\n");
  print_mat(B[1], nmo, nmo, outfile);
  fprintf(outfile, "\tMO-basis MuZ Integrals:\n");
  print_mat(B[2], nmo, nmo, outfile);
  */

  /* Sort the B's into vector storage */
  B0 = block_matrix(3, num_ai);
  for(coord=0; coord < 3; coord++) {
    for(asym=0,AI=0; asym < nirreps; asym++) {

      afirst = vfirst[asym];
      alast = vlast[asym];

      for(a=afirst; a <= alast; a++) {

	for(isym = 0; isym < nirreps; isym++) {
	  ifirst = ofirst[isym];
	  ilast = olast[isym];

	  for(i=ifirst; i <= ilast; i++,AI++) {

	    B0[coord][AI] = B[coord][a][i];
	  }
	}
      }
    }
  }

  ipiv = init_int_array(num_ai);

  /* Solve the CPHF equations */
  Acopy = block_matrix(num_ai, num_ai); /* keep a copy of A */
  memcpy(Acopy[0], A[0], num_ai*num_ai*sizeof(double));

  for(coord=0; coord < 3; coord++) {
    error = C_DGESV(num_ai, 1, &(A[0][0]), num_ai, &(ipiv[0]), &(B0[coord][0]), num_ai);

    /* Recopy A because DGESV corrupts it */
    memcpy(A[0], Acopy[0], num_ai*num_ai*sizeof(double));
  }

  /* Sort the U matrices to matrix form */
  for(coord=0; coord < 3; coord++) {
    for(asym=0,AI=0; asym < nirreps; asym++) {

      afirst = vfirst[asym];
      alast = vlast[asym];

      for(a=afirst; a <= alast; a++) {

	for(isym=0; isym < nirreps; isym++) {

	  ifirst = ofirst[isym];
	  ilast = olast[isym];

	  for(i=ifirst; i <= ilast; i++,AI++) {

	    U[coord][a][i] = B0[coord][AI];
	  }
	}
      }
    }
  }

  /* Dump the U matrices out to disk */
  label = (char *) malloc(PSIO_KEYLEN * sizeof(char));
  psio_open(PSIF_CPHF, 1);
  for(coord=0; coord < 3; coord++) {
    sprintf(label, "UF(%d)", coord);
    psio_write_entry(PSIF_CPHF, label, (char *) &(U[coord][0][0]), nmo*nmo*sizeof(double));
    for(i=0; i < PSIO_KEYLEN; i++) label[i] = '\0';

    /*
      fprintf(outfile, "\nUF[%d] Matrix (MO):\n", coord);
      print_mat(U[coord], nmo, nmo, outfile);
    */
  }
  psio_close(PSIF_CPHF, 1);

  for(coord=0; coord < 3; coord++)
    free_block(B[coord]);

  free(B);
  free_block(B0); 

  free_block(Acopy);
}
