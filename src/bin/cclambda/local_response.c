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
  int i, j, ij, a, k, max, complete, *boolean, *rank;
  int nao, nso, nmo,noei_ao;
  double **TMP, *scratch, **X;
  double **MU, **Z;
  double **C, **usotao;
  double polar, polar_i, polar_k, polar_k_check, value, *polar_mo, *polar_mo_check, **polar_atom;
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
  polar_mo_check = init_array(moinfo.occpi[0]);
  polar_mo = init_array(moinfo.occpi[0]);
  polar_atom = block_matrix(moinfo.occpi[0],natom);
  boolean = init_int_array(natom);
  rank = (int *) malloc(natom * sizeof(int *));

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
    polar_mo_check[i] = 0.0;
    polar_mo[i] = 0.0;
    for(k=0; k < natom; k++) {
      polar_atom[i][k] = 0.0;
      for(a=aostart[k]; a <= aostop[k]; a++) {
	value = 4.0 * Z[a][i] * MU[a][i];
	polar += value;
	polar_mo[i] += value;
	polar_mo_check[i] += fabs(value);
	polar_atom[i][k] += fabs(value);
      }
    }
  }
  /*   fprintf(outfile, "\n\tTotal alpha_xx = %20.12f\n", polar); */

  for(i=0; i < moinfo.occpi[0]; i++) {
    /* Rank the atomic contributions to the orbital's polarizability */
    for(j=0; j < natom; j++) {
      rank[j] = 0;
      boolean[j] = 0;
    }
    for(j=0,max=0; j < natom; j++) /* find the overall maximum */
      if(fabs(polar_atom[i][j]) >= fabs(polar_atom[i][max])) max = j;
    rank[0] = max;
    boolean[max] = 1;
    for(j=1; j < natom; j++) {
      max = 0;
      while(boolean[max]) max++; /* find an unused max */
      for(k=0; k < natom; k++)
	if((fabs(polar_atom[i][k]) >= fabs(polar_atom[i][max])) && !boolean[k]) max = k;
      rank[j] = max; 
      boolean[max] = 1;
    }

    /*     for(k=0; k < natom; k++) */
    /*       fprintf(outfile, "\tAlpha_xx_check[%d][%d] = %15.12lf\n", i, rank[k], polar_atom[i][rank[k]]); */

    /* Response domains for Alpha_xx */
    polar_i = 0.0;
    complete = 0;
    /*     fprintf(outfile, "\talpha_xx_mo[%d] = %20.12f\n", i, polar_mo[i]); */
    /*     fprintf(outfile, "\talpha_xx_mo_check[%d] = %20.12f\n", i, polar_mo_check[i]); */
    for(k=0; k < natom; k++) {
      polar_k = 0.0;
      polar_k_check = 0.0;
      for(a=aostart[rank[k]]; a <= aostop[rank[k]]; a++) {
	value = Z[a][i] * MU[a][i];
	polar_k += value;
	polar_k_check += fabs(value);
      }
      polar_k *= 4.0;
      polar_k_check *= 4.0;
      polar_i += polar_k_check;
      /*       fprintf(outfile, "\talpha_xx[%d][%d] = %20.12f.", i, rank[k], polar_k); */
      if (!complete) {
	if(!domain[i][rank[k]]) {
	  /* 	  fprintf(outfile, "  Adding [%d][%d] to domain list.", i, rank[k]); */
	  domain[i][rank[k]] = 1;
	  domain_len[i]++;
	}
	/* 	else if(domain[i][rank[k]]) */
	/* 	  fprintf(outfile, "  [%d][%d] is already in domain list.", i, rank[k]); */
	if(fabs((polar_mo_check[i]-polar_i)/polar_mo_check[i]) <= local.cphf_cutoff) {
	  complete = 1;
	  if(fabs(polar_atom[i][k]-polar_atom[i][k+1] <= 0.0001)) {
	    if(!domain[i][rank[k]]) {
	      /* 	      fprintf(outfile, "  Adding [%d][%d] to domain list.", i, rank[k]); */
	      domain[i][rank[k]] = 1;
	      domain_len[i]++;
	    }
	  }
	}
      }
      /*       else if(domain[i][rank[k]]) { */
      /* 	fprintf(outfile, "  [%d][%d] is already in domain list.", i, rank[k]); */
      /*       } */
      /*       fprintf(outfile, " \n"); */
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
    polar_mo_check[i] = 0.0;
    polar_mo[i] = 0.0;
    for(k=0; k < natom; k++) {
      polar_atom[i][k] = 0.0;
      for(a=aostart[k]; a <= aostop[k]; a++) {
	value = 4.0 * Z[a][i] * MU[a][i];
	polar += value;
	polar_mo[i] += value;
	polar_mo_check[i] += fabs(value);
	polar_atom[i][k] += fabs(value);
      }
    }
  }
  /*   fprintf(outfile, "\n\tTotal alpha_yy = %20.12f\n", polar); */

  for(i=0; i < moinfo.occpi[0]; i++) {
    /* Rank the atomic contributions to the orbital's polarizability */
    for(j=0; j < natom; j++) {
      rank[j] = 0;
      boolean[j] = 0;
    }
    for(j=0,max=0; j < natom; j++) /* find the overall maximum */
      if(fabs(polar_atom[i][j]) >= fabs(polar_atom[i][max])) max = j;
    rank[0] = max; boolean[max] = 1;
    for(j=1; j < natom; j++) {
      max = 0;
      while(boolean[max]) max++; /* find an unused max */
      for(k=0; k < natom; k++)
	if((fabs(polar_atom[i][k]) >= fabs(polar_atom[i][max])) && !boolean[k]) max = k;
      rank[j] = max; 
      boolean[max] = 1;
    }

    /*     for(k=0; k < natom; k++) */
    /*       fprintf(outfile, "\tAlpha_yy_check[%d][%d] = %15.12lf\n", i, rank[k], polar_atom[i][rank[k]]); */

    /* Response domains for Alpha_yy */
    polar_i = 0.0;
    complete = 0;
    /*     fprintf(outfile, "\talpha_yy_mo[%d] = %20.12f\n", i, polar_mo[i]); */
    /*     fprintf(outfile, "\talpha_yy_mo_check[%d] = %20.12f\n", i, polar_mo_check[i]); */
    for(k=0; k < natom; k++) {
      polar_k = 0.0;
      polar_k_check = 0.0;
      for(a=aostart[rank[k]]; a <= aostop[rank[k]]; a++) {
	value = Z[a][i] * MU[a][i];
	polar_k += value;
	polar_k_check += fabs(value);
      }
      polar_k *= 4.0;
      polar_k_check *= 4.0;
      polar_i += polar_k_check;
      /*       fprintf(outfile, "\talpha_yy[%d][%d] = %20.12f.", i, rank[k], polar_k); */
      if (!complete) {
	if(!domain[i][rank[k]]) {
	  /* 	  fprintf(outfile, "  Adding [%d][%d] to domain list.", i, rank[k]); */
	  domain[i][rank[k]] = 1;
	  domain_len[i]++;
	}
	/* 	else if(domain[i][rank[k]]) { */
	/* 	  fprintf(outfile, "  [%d][%d] is already in domain list.", i, rank[k]); */
	/* 	} */
	if(fabs((polar_mo_check[i]-polar_i)/polar_mo_check[i]) <= local.cphf_cutoff) {
	  complete = 1;
	  if(fabs(polar_atom[i][k]-polar_atom[i][k+1] <= 0.0001)) {
	    if(!domain[i][rank[k]]) {
	      /* 	      fprintf(outfile, "  Adding [%d][%d] to domain list.", i, rank[k]); */
	      domain[i][rank[k]] = 1;
	      domain_len[i]++;
	    }
	  }
	}
      }
      /*       else if(domain[i][rank[k]]) { */
      /* 	fprintf(outfile, "  [%d][%d] is already in domain list.", i, rank[k]); */
      /*       } */
      /*       fprintf(outfile, " \n"); */
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
    polar_mo[i] = 0.0;
    polar_mo_check[i] = 0.0;
    for(k=0; k < natom; k++) {
      polar_atom[i][k] = 0.0;
      for(a=aostart[k]; a <= aostop[k]; a++) {
	value = 4.0 * Z[a][i] * MU[a][i];
	polar += value;
	polar_mo[i] += value;
	polar_mo_check[i] += fabs(value);
	polar_atom[i][k] += fabs(value);
      }
    }
  }
  /*   fprintf(outfile, "\n\tTotal alpha_zz = %20.12f\n", polar); */

  for(i=0; i < moinfo.occpi[0]; i++) {
    /* Rank the atomic contributions to the orbital's polarizability */
    for(j=0; j < natom; j++) {
      rank[j] = 0;
      boolean[j] = 0;
    }
    for(j=0,max=0; j < natom; j++) /* find the overall maximum */
      if(fabs(polar_atom[i][j]) >= fabs(polar_atom[i][max])) max = j;
    rank[0] = max; boolean[max] = 1;
    for(j=1; j < natom; j++) {
      max = 0;
      while(boolean[max]) max++; /* find an unused max */
      for(k=0; k < natom; k++)
	if((fabs(polar_atom[i][k]) >= fabs(polar_atom[i][max])) && !boolean[k]) max = k;
      rank[j] = max; 
      boolean[max] = 1;
    }

    /*      for(k=0; k < natom; k++) */
    /*       fprintf(outfile, "\tAlpha_zz_check[%d][%d] = %15.12lf\n", i, rank[k], polar_atom[i][rank[k]]); */

    /* Response domains for Alpha_zz */
    polar_i = 0.0;
    complete = 0;
    /*     fprintf(outfile, "\talpha_zz_mo[%d] = %20.12f\n", i, polar_mo[i]); */
    /*     fprintf(outfile, "\talpha_zz_mo_check[%d] = %20.12f\n", i, polar_mo_check[i]); */
    for(k=0; k < natom; k++) {
      polar_k = 0.0;
      polar_k_check = 0.0;
      for(a=aostart[rank[k]]; a <= aostop[rank[k]]; a++) {
	value = Z[a][i] * MU[a][i];
	polar_k += value;
	polar_k_check += fabs(value);
      }
      polar_k *= 4.0;
      polar_k_check *= 4.0;
      polar_i += polar_k_check;
      /*       fprintf(outfile, "\talpha_zz[%d][%d] = %20.12f.", i, rank[k], polar_k); */
      if (!complete) {
	if(!domain[i][rank[k]]) {
	  fprintf(outfile, "  Adding [%d][%d] to domain list.", i, rank[k]);
	  domain[i][rank[k]] = 1;
	  domain_len[i]++;
	}
	/* 	else if(domain[i][rank[k]]) { */
	/* 	  fprintf(outfile, "  [%d][%d] is already in domain list.", i, rank[k]); */
	/* 	} */
	if(fabs((polar_mo_check[i]-polar_i)/polar_mo_check[i]) <= local.cphf_cutoff) {
	  complete = 1;
	  if(fabs(polar_atom[i][k]-polar_atom[i][k+1] <= 0.0001)) {
	    if(!domain[i][rank[k]]) {
	      /* 	      fprintf(outfile, "  Adding [%d][%d] to domain list.", i, rank[k]); */
	      domain[i][rank[k]] = 1;
	      domain_len[i]++;
	    }
	  }
	}
      }
      /*       else if(domain[i][rank[k]]) { */
      /* 	fprintf(outfile, "  [%d][%d] is already in domain list.", i, rank[k]); */
      /*       } */
      /*       fprintf(outfile, " \n"); */
    }
  }

  dpd_file2_mat_close(&U);
  dpd_file2_close(&U);

  free_block(TMP);
  free_block(MU);

  free_block(polar_atom);
  free_block(C);

  free(boolean);
  free(polar_mo_check);
  free(polar_mo);
  free(rank);
}
