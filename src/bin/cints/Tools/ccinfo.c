#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<libciomr/libciomr.h>
#include<libpsio/psio.h>
#include<libdpd/dpd.h>
#include<libchkpt/chkpt.h>
#include<libint/libint.h>
#include<libqt/qt.h>
#include <psifiles.h>
#include <ccfiles.h>

#include"moinfo.h"
#include"moinfo_corr.h"

#include"defines.h"
#define EXTERN
#include"global.h"
#include"small_fns.h"

static int *cachefiles, **cachelist;

void init_ccinfo()
{
  int h, ij, ab, row, col, row_offset, col_offset;
  int nactive, nocc, nvirt;
  int *occpi, *occ_sym;
  int *virtpi, *vir_sym;
  double **A;
  double **T2_MO, **T2_AO, **T2;
  dpdbuf4 tau;

  /* grab some basic data from CC_INFO */
  psio_open(CC_INFO, PSIO_OPEN_OLD);
  psio_read_entry(CC_INFO, "No. of Active Orbitals", (char *) &(nactive),
                  sizeof(int));
  occpi = init_int_array(Symmetry.nirreps);
  virtpi = init_int_array(Symmetry.nirreps);
  occ_sym = init_int_array(nactive);
  vir_sym = init_int_array(nactive);
  psio_read_entry(CC_INFO, "Active Occ Orbs Per Irrep",
                  (char *) occpi, sizeof(int)*Symmetry.nirreps);
  psio_read_entry(CC_INFO, "Active Virt Orbs Per Irrep",
                  (char *) virtpi, sizeof(int)*Symmetry.nirreps);
  psio_read_entry(CC_INFO, "Active Occ Orb Symmetry",
                  (char *) occ_sym, sizeof(int)*Symmetry.nirreps);
  psio_read_entry(CC_INFO, "Active Virt Orbs Symmetry",
                  (char *) vir_sym, sizeof(int)*Symmetry.nirreps);
  psio_close(CC_INFO, 1);

  nocc = 0; nvirt = 0;
  for(h=0; h < Symmetry.nirreps; h++) {
    nocc += occpi[h];  nvirt += virtpi[h];
  }

  cachefiles = init_int_array(PSIO_MAXUNIT);
  /* assuming no caching for DPD here */
  cachelist = init_int_matrix(12,12);

  init_moinfo();
  init_moinfo_corr();

  /*--- Initialize DPD library ---*/
  dpd_init(1, Symmetry.nirreps, UserOptions.memory, 0, cachefiles, cachelist, NULL, 
	   2, occpi, occ_sym, virtpi, vir_sym);

  /* Grab the MO-basis T2's provided by ccenergy */
  T2_MO = block_matrix(nocc*nocc,BasisSet.num_ao*BasisSet.num_ao);
  dpd_buf4_init(&tau, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
  for(h=0,row=0,col=0; h < Symmetry.nirreps; h++) {
    dpd_buf4_mat_irrep_init(&tau, h);
    dpd_buf4_mat_irrep_rd(&tau, h);

    if(h) { 
      row_offset += tau.params->rowtot[h-1];  
      col_offset += tau.params->coltot[h-1];
    }

    for(row=0; row < tau.params->rowtot[h]; row++) {
      for(col=0; col < tau.params->coltot[h]; col++) {
        T2_MO[row_offset+row][col_offset+col] = tau.matrix[h][row][col];
      }
    }

  dpd_buf4_mat_irrep_close(&tau, h);
  }
  dpd_buf4_close(&tau);

  /* half-transform the T2's */
  A = block_matrix(nvirt, BasisSet.num_ao);
  T2_AO = T2_MO;
  for(ij=0; ij < nocc*nocc; ij++) {
    
    C_DGEMM('n','n',nvirt,BasisSet.num_ao,nvirt,1.0,
            T2_MO[ij],BasisSet.num_ao*BasisSet.num_ao,
            MOInfo.scf_evec_uocc[0][0],BasisSet.num_ao,0.0,A[0],BasisSet.num_ao);
    C_DGEMM('t','n',BasisSet.num_ao,BasisSet.num_ao,nvirt,1.0,
            MOInfo.scf_evec_uocc[0][0],BasisSet.num_ao,A[0],BasisSet.num_ao,
            0.0,T2_AO[ij],BasisSet.num_ao*BasisSet.num_ao);
  }
  free_block(A);

  /* transpose T2 to AO*AO x MO*MO */
  T2 = block_matrix(BasisSet.num_ao*BasisSet.num_ao,nocc*nocc);
  for(ij=0; ij < nocc*nocc; ij++)
    for(ab=0; ab < BasisSet.num_ao*BasisSet.num_ao; ab++) 
      T2[ab][ij] = T2_AO[ij][ab];
  free_block(T2_AO);

  CCInfo.T2_s = T2;
  CCInfo.T2_t = block_matrix(BasisSet.num_ao*BasisSet.num_ao,nocc*nocc);
  CCInfo.nocc = nocc;
  CCInfo.nvirt = nvirt;

  free(occpi);  free(occ_sym);
  free(virtpi); free(vir_sym);
      
  return;
}


void cleanup_ccinfo()
{
  free_block(CCInfo.T2_s);
  free_block(CCInfo.T2_t);

  free_int_matrix(cachelist,12);
  free(cachefiles);
  dpd_close(1);

  return;
}


