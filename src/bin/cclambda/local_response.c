#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libiwl/iwl.h>
#include <libchkpt/chkpt.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

void local_response(int **domain, int *domain_len, int natom, int *aostart, int *aostop)
{
  int i, j, ij, a, k;
  int nao, nso, nmo,noei_ao;
  double **TMP, *scratch, **X;
  double **MU, **Z;
  double **C, **usotao;
  double polar, polar_i, polar_k;
  dpdfile2 U;
  psio_address next;
  int stat;

  nao = moinfo.nao;
  nso = moinfo.nso;
  nmo = moinfo.nmo;
  noei_ao = nao*(nao+1)/2;

  /* grab the occupied MOs */
  next = PSIO_ZERO;
  C = block_matrix(moinfo.nso, moinfo.occpi[0]);
  psio_read(CC_INFO, "RHF/ROHF Active Occupied Orbitals", (char *) C[0],
	    nso*moinfo.occpi[0]*sizeof(double), next, &next);

  /* grab the usotao matrix */
  chkpt_init(PSIO_OPEN_OLD);
  usotao = chkpt_rd_usotao();
  chkpt_close();

  MU = block_matrix(nso,moinfo.occpi[0]);
  TMP = block_matrix(nao,nao);
  X = block_matrix(nao, nao);
  scratch = init_array(noei_ao);
  Z = block_matrix(nso, moinfo.occpi[0]);

  stat = iwl_rdone(PSIF_OEI, PSIF_AO_MX, scratch, noei_ao, 0, 0, outfile);
  for(i=0, ij=0; i < nao; i++)
    for(j=0; j <= i; j++,ij++)
      TMP[i][j] = TMP[j][i] = scratch[ij];

  C_DGEMM('n','t',nao,nso,nao,1,&(TMP[0][0]),nao,&(usotao[0][0]),nao,
	  0,&(X[0][0]),nao);
  C_DGEMM('n','n',nso,nso,nao,1,&(usotao[0][0]),nao,&(X[0][0]),nao,
	  0,&(TMP[0][0]),nao);

  C_DGEMM('n', 'n', nso, moinfo.occpi[0], nso, 1, &(TMP[0][0]), nao, &(C[0][0]), moinfo.occpi[0],
	  0, &(MU[0][0]), moinfo.occpi[0]);

  dpd_file2_init(&U, CC_OEI, 0, 1, 0, "CPHF U_X_AI");
  dpd_file2_mat_init(&U);
  dpd_file2_mat_rd(&U);

  C_DGEMM('n', 'n', nso,moinfo.occpi[0],moinfo.virtpi[0], 1, &(moinfo.C[0][0][0]), moinfo.virtpi[0],
	  &(U.matrix[0][0][0]), moinfo.occpi[0], 0, &(Z[0][0]), moinfo.occpi[0]);

  polar = 0.0;
  for(i=0; i < moinfo.occpi[0]; i++) {
    for(k=0; k < natom; k++) {
      for(a=aostart[k]; a <= aostop[k]; a++) {
	polar += 4.0 * Z[a][i] * MU[a][i];
      }
    }
  }
  fprintf(outfile, "\tTotal alpha_xx = %20.12f\n", polar);

  for(i=0; i < moinfo.occpi[0]; i++) {
    for(k=0; k < natom; k++) {
      polar_k = 0.0;
      for(a=aostart[k]; a <= aostop[k]; a++) {
	polar_k += Z[a][i] * MU[a][i];
      }
      polar_k *= 4.0;
      fprintf(outfile, "\talpha_xx[%d][%d] = %20.12f.", i, k, polar_k);
      if(fabs(polar_k/polar) > local.cphf_cutoff && !domain[i][k]) {
	fprintf(outfile, "  Adding [%d][%d] to domain list.\n", i, k);
	domain[i][k] = 1;
	domain_len[i]++;
      }
      else if(fabs(polar_k/polar) > local.cphf_cutoff && domain[i][k]) {
	fprintf(outfile, "  [%d][%d] is already in domain list.\n", i, k);
      }
      else fprintf(outfile, "\n");
    }
  }

  dpd_file2_mat_close(&U);
  dpd_file2_close(&U);

  stat = iwl_rdone(PSIF_OEI, PSIF_AO_MY, scratch, noei_ao, 0, 0, outfile);
  for(i=0, ij=0; i < nao; i++)
    for(j=0; j <= i; j++,ij++)
      TMP[i][j] = TMP[j][i] = scratch[ij];

  C_DGEMM('n','t',nao,nso,nao,1,&(TMP[0][0]),nao,&(usotao[0][0]),nao,
	  0,&(X[0][0]),nao);
  C_DGEMM('n','n',nso,nso,nao,1,&(usotao[0][0]),nao,&(X[0][0]),nao,
	  0,&(TMP[0][0]),nao);

  C_DGEMM('n', 'n', nso, moinfo.occpi[0], nso, 1, &(TMP[0][0]), nao, &(C[0][0]), moinfo.occpi[0],
	  0, &(MU[0][0]), moinfo.occpi[0]);

  dpd_file2_init(&U, CC_OEI, 0, 1, 0, "CPHF U_Y_AI");
  dpd_file2_mat_init(&U);
  dpd_file2_mat_rd(&U);

  C_DGEMM('n', 'n', nso,moinfo.occpi[0],moinfo.virtpi[0], 1, &(moinfo.C[0][0][0]), moinfo.virtpi[0],
	  &(U.matrix[0][0][0]), moinfo.occpi[0], 0, &(Z[0][0]), moinfo.occpi[0]);

  polar = 0.0;
  for(i=0; i < moinfo.occpi[0]; i++) {
    for(k=0; k < natom; k++) {
      for(a=aostart[k]; a <= aostop[k]; a++) {
	polar += 4.0 * Z[a][i] * MU[a][i];
      }
    }
  }
  fprintf(outfile, "\tTotal alpha_yy = %20.12f\n", polar);

  for(i=0; i < moinfo.occpi[0]; i++) {
    for(k=0; k < natom; k++) {
      polar_k = 0.0;
      for(a=aostart[k]; a <= aostop[k]; a++) {
	polar_k += Z[a][i] * MU[a][i];
      }
      polar_k *= 4.0;
      fprintf(outfile, "\talpha_yy[%d][%d] = %20.12f.", i, k, polar_k);
      if(fabs(polar_k/polar) > local.cphf_cutoff && !domain[i][k]) {
	fprintf(outfile, "   Adding [%d][%d] to domain list.\n", i, k);
	domain[i][k] = 1;
	domain_len[i]++;
      }
      else if(fabs(polar_k/polar) > local.cphf_cutoff && domain[i][k]) {
	fprintf(outfile, "  [%d][%d] is already in domain list.\n", i, k);
      }
      else fprintf(outfile, "\n");
    }
  }

  dpd_file2_mat_close(&U);
  dpd_file2_close(&U);

  stat = iwl_rdone(PSIF_OEI, PSIF_AO_MZ, scratch, noei_ao, 0, 0, outfile);
  for(i=0, ij=0; i < nao; i++)
    for(j=0; j <= i; j++,ij++)
      TMP[i][j] = TMP[j][i] = scratch[ij];

  C_DGEMM('n','t',nao,nso,nao,1,&(TMP[0][0]),nao,&(usotao[0][0]),nao,
	  0,&(X[0][0]),nao);
  C_DGEMM('n','n',nso,nso,nao,1,&(usotao[0][0]),nao,&(X[0][0]),nao,
	  0,&(TMP[0][0]),nao);

  C_DGEMM('n', 'n', nso, moinfo.occpi[0], nso, 1, &(TMP[0][0]), nao, &(C[0][0]), moinfo.occpi[0],
	  0, &(MU[0][0]), moinfo.occpi[0]);

  dpd_file2_init(&U, CC_OEI, 0, 1, 0, "CPHF U_Z_AI");
  dpd_file2_mat_init(&U);
  dpd_file2_mat_rd(&U);

  C_DGEMM('n', 'n', nso,moinfo.occpi[0],moinfo.virtpi[0], 1, &(moinfo.C[0][0][0]), moinfo.virtpi[0],
	  &(U.matrix[0][0][0]), moinfo.occpi[0], 0, &(Z[0][0]), moinfo.occpi[0]);

  polar = 0.0;
  for(i=0; i < moinfo.occpi[0]; i++) {
    for(k=0; k < natom; k++) {
      for(a=aostart[k]; a <= aostop[k]; a++) {
	polar += 4.0 * Z[a][i] * MU[a][i];
      }
    }
  }
  fprintf(outfile, "\tTotal alpha_zz = %20.12f\n", polar);

  for(i=0; i < moinfo.occpi[0]; i++) {
    for(k=0; k < natom; k++) {
      polar_k = 0.0;
      for(a=aostart[k]; a <= aostop[k]; a++) {
	polar_k += Z[a][i] * MU[a][i];
      }
      polar_k *= 4.0;
      fprintf(outfile, "\talpha_zz[%d][%d] = %20.12f.", i, k, polar_k);
      if(fabs(polar_k/polar) > local.cphf_cutoff && !domain[i][k]) {
	fprintf(outfile, "   Adding [%d][%d] to domain list.\n", i, k);
	domain[i][k] = 1;
	domain_len[i]++;
      }
      else if(fabs(polar_k/polar) > local.cphf_cutoff && domain[i][k]) {
	fprintf(outfile, "   [%d][%d] is already in domain list.\n", i, k);
      }
      else fprintf(outfile, "\n");
    }
  }

  dpd_file2_mat_close(&U);
  dpd_file2_close(&U);

  free_block(TMP);
  free_block(MU);

  free_block(C);
}
