#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#include <physconst.h>
#include "globals.h"

void init_io(int argc, char *argv[]);
void title(void);
void exit_io(void);
int **cacheprep_uhf(int level, int *cachefiles);
int **cacheprep_rhf(int level, int *cachefiles);
void cachedone_rhf(int **cachelist);
void cachedone_uhf(int **cachelist);
void get_moinfo(void);
void get_params(void);
void cleanup(void);
void build_A(void);
void d_corr(void);
void local_init(void);
void local_done(void);
void amp_write_T1(dpdfile2 *T1, int length, FILE *outfile);

int main(int argc, char *argv[])
{
  char lbl[32];
  int **cachelist, *cachefiles, *ao2atom;
  int h, i, j, jj, k, nroot, a, max_j, max_a;
  double *evals, *dcorr, *weakp, value, d_value, **B_PAO;
  double *Bt, max_val, test_val, sum_val;
  int *spin, *symm, count, ivalue;
  dpdfile2 B;

  init_io(argc, argv);
  title();

  get_moinfo();
  get_params();

  cachefiles = init_int_array(PSIO_MAXUNIT);

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/
    cachelist = cacheprep_rhf(0, cachefiles);

    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, 
        cachelist, NULL, 2, moinfo.occpi, moinfo.occ_sym, 
        moinfo.virtpi, moinfo.vir_sym);
  }
  else if(params.ref == 2) { /** UHF **/
    cachelist = cacheprep_uhf(0, cachefiles);

    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, 
        cachelist, NULL, 4, moinfo.aoccpi, moinfo.aocc_sym, moinfo.avirtpi,
        moinfo.avir_sym, moinfo.boccpi, moinfo.bocc_sym, moinfo.bvirtpi, moinfo.bvir_sym);
  }

  /* build components of CIS response matrix, A */
  build_A();

  /* diagonalize A within each irrep */
  diag();

  /* moved up with impunity??? */
  if(params.local) local_init();

  if(params.local) {
    ao2atom = init_int_array(local.nso);
    for(i=0; i < local.natom; i++)
      for(j=local.aostart[i]; j <= local.aostop[i]; j++) ao2atom[j] = i;

  /*  for(i=0; i < local.nso; i++)
      fprintf(outfile,"ao2atom[%d] = %d\n", i, ao2atom[i]);
      */
  }

  /* print eigenvalues and largest components of eigenvectors in ascending order */
  if(params.ref == 0) {

    fprintf(outfile, "\tRHF-CIS Singlet Excitation Energies:\n");
    fprintf(outfile, "\t------------------------------------\n\n");
    fprintf(outfile, "\tRoot Irrep       Hartree          eV          cm-1  \n");
    fprintf(outfile, "\t---- ----- ------------------  ---------  ----------\n");
    for(h=0; h < moinfo.nirreps; h++) {
      for(i=0; i < params.rpi[h]; i++) {
        value = moinfo.singlet_evals[h][i];

        fprintf(outfile,   "\t%4d   %3s %18.14f  %9.5f  %10.2f\n", i, moinfo.labels[h],
            value, value*_hartree2ev, value*_hartree2wavenumbers);

        fprintf(outfile, "\nLargest components of singlet excited wave function #%d/#%d:\n",
            h, i);
        sprintf(lbl, "BIA(%d)[%d] singlet", i, h);
        dpd_file2_init(&B, CC_OEI, h, 0, 1, lbl);
        amp_write_T1(&B, 5, outfile);
        dpd_file2_close(&B);
        
        fprintf(outfile,"\nLargest components in projected AO virtual basis\n");
        if (params.local) {
          /* Transform the virtuals to the redundant projected virtual basis */
          sprintf(lbl, "BIA(%d)[%d] singlet", i, h);
          dpd_file2_init(&B, CC_OEI, h, 0, 1, lbl);
          dpd_file2_mat_init(&B);
          dpd_file2_mat_rd(&B);
          B_PAO = block_matrix(local.nocc, local.nso);
          C_DGEMM('n','n', local.nocc, local.nso, local.nvir, 1.0, &(B.matrix[0][0][0]),
              local.nvir, &(local.U[0][0]), local.nso, 0.0, &(B_PAO[0][0]), local.nso);
          dpd_file2_mat_close(&B);
          dpd_file2_close(&B);
  
          // print_mat(B_PAO, local.nocc, local.nso, outfile);
          fprintf(outfile,"\t       occ  vir  atom        amplitude\n");
          sum_val = 0.0;
          for (jj=0; jj< local.nocc*local.nso; ++jj) {
            max_val = 0.0;
            for (j=0; j<local.nocc; ++j) {
              for (a=0; a<local.nso; ++a) {
                test_val = fabs(B_PAO[j][a]);
                if (test_val > max_val) {
	                max_j = j;
	                max_a = a;
	                max_val = test_val;
                }
              }
            }
            if (sum_val < local.amp_print_cutoff) {
              fprintf(outfile,"\t      %3d  %3d %4d %20.10f\n",
                max_j, max_a, ao2atom[max_a], B_PAO[max_j][max_a]);
              sum_val += fabs(B_PAO[max_j][max_a]);
              B_PAO[max_j][max_a] = 0.0;
            }
	        }
          free_block(B_PAO);
          fprintf(outfile, "\n");
        }
      }
    }
    fprintf(outfile, "\n");

    fprintf(outfile, "\tRHF-CIS Triplet Excitation Energies:\n");
    fprintf(outfile, "\t------------------------------------\n\n");
    fprintf(outfile, "\tRoot Irrep       Hartree          eV          cm-1  \n");
    fprintf(outfile, "\t---- ----- ------------------  ---------  ----------\n");
    for(h=0; h < moinfo.nirreps; h++) {
      for(i=0; i < params.rpi[h]; i++) {
        value = moinfo.triplet_evals[h][i];

        fprintf(outfile,   "\t%4d   %3s %18.14f  %9.5f  %10.2f\n", i, moinfo.labels[h],
            value, value*_hartree2ev, value*_hartree2wavenumbers);

        fprintf(outfile, "\nLargest components of triplet excited wave function #%d/#%d:\n", h, i);
        sprintf(lbl, "BIA(%d)[%d] triplet", i, h);
        dpd_file2_init(&B, CC_OEI, h, 0, 1, lbl);
        amp_write_T1(&B, 5, outfile);
        dpd_file2_close(&B);
        fprintf(outfile, "\n");

        if (params.local) {
          /* Transform the virtuals to the redundant projected virtual basis */
          sprintf(lbl, "BIA(%d)[%d] triplet", i, h);
          dpd_file2_init(&B, CC_OEI, h, 0, 1, lbl);
          dpd_file2_mat_init(&B);
          dpd_file2_mat_rd(&B);
          B_PAO = block_matrix(local.nocc, local.nso);
          C_DGEMM('n','n', local.nocc, local.nso, local.nvir, 1.0, &(B.matrix[0][0][0]),
              local.nvir, &(local.U[0][0]), local.nso, 0.0, &(B_PAO[0][0]), local.nso);
          dpd_file2_mat_close(&B);
          dpd_file2_close(&B);
  
          // print_mat(B_PAO, local.nocc, local.nso, outfile);
          fprintf(outfile,"\t       occ  vir  atom        amplitude\n");
          sum_val = 0.0;
          for (jj=0; jj< local.nocc*local.nso; ++jj) {
            max_val = 0.0;
            for (j=0; j<local.nocc; ++j) {
              for (a=0; a<local.nso; ++a) {
                test_val = fabs(B_PAO[j][a]);
                if (test_val > max_val) {
                  max_j = j;
                  max_a = a;
                  max_val = test_val;
                }
              }
            }
            if (sum_val < local.amp_print_cutoff) {
              fprintf(outfile,"\t      %3d  %3d %4d %20.10f\n",
                max_j, max_a, ao2atom[max_a], B_PAO[max_j][max_a]);
              sum_val += fabs(B_PAO[max_j][max_a]);
              B_PAO[max_j][max_a] = 0.0;
            }
          }
          free_block(B_PAO);
          fprintf(outfile, "\n");
        }
      }
    }
    fprintf(outfile, "\n");

  }
  free(ao2atom);

  fflush(outfile);

  /* compute the (D) correction to each CIS singlet excitation energy */
  // move up with impunity???  if(params.local) local_init();

  d_corr();

  fprintf(outfile, "\n");

  if(params.local) local_done();

  /* print corrected eigenvalues in ascending order */
  if(params.ref == 0) {

    fprintf(outfile, "\tRHF-CIS(D) Singlet Corrections:\n");
    fprintf(outfile, "\t-------------------------------\n");
    fprintf(outfile, "\tRoot Irrep       Hartree          eV          cm-1  \n");
    fprintf(outfile, "\t---- ----- ------------------  ---------  ----------\n");
    for(h=0; h < moinfo.nirreps; h++) {
      for(i=0; i < params.rpi[h]; i++) {
        value = moinfo.singlet_d[h][i];

        fprintf(outfile,   "\t%4d   %3s %18.14f  %9.5f  %10.2f\n", i, moinfo.labels[h],
            value, value*_hartree2ev, value*_hartree2wavenumbers);
      }
    }
    fprintf(outfile, "\n");

    fprintf(outfile, "\tRHF-CIS(D) Singlet Excitation Energies:\n");
    fprintf(outfile, "\t---------------------------------------\n\n");
    fprintf(outfile, "\tRoot Irrep       Hartree          eV          cm-1  \n");
    fprintf(outfile, "\t---- ----- ------------------  ---------  ----------\n");
    for(h=0; h < moinfo.nirreps; h++) {
      for(i=0; i < params.rpi[h]; i++) {
        value = moinfo.singlet_evals[h][i];
        value += moinfo.singlet_d[h][i];

        fprintf(outfile,   "\t%4d   %3s %18.14f  %9.5f  %10.2f\n", i, moinfo.labels[h],
            value, value*_hartree2ev, value*_hartree2wavenumbers);
      }
    }
    fprintf(outfile, "\n");

    if(params.local) {
      fprintf(outfile, "\tRHF-CIS(D) Weak Pair Corrections:\n");
      fprintf(outfile, "\t---------------------------------\n");
      fprintf(outfile, "\tRoot Irrep       Hartree          eV          cm-1  \n");
      fprintf(outfile, "\t---- ----- ------------------  ---------  ----------\n");
      for(h=0; h < moinfo.nirreps; h++) {
        for(i=0; i < params.rpi[h]; i++) {
          value = moinfo.singlet_weakp[h][i];

          fprintf(outfile,   "\t%4d   %3s %18.14f  %9.5f  %10.2f\n", i, moinfo.labels[h],
              value, value*_hartree2ev, value*_hartree2wavenumbers);
        }
      }
      fprintf(outfile, "\n");
    }
  }

  else if(params.ref == 2) {

    fprintf(outfile, "\tUHF-CIS Excitation Energies:\n");
    fprintf(outfile, "\t----------------------------\n\n");
    fprintf(outfile, "\tRoot Irrep       Hartree          eV          cm-1\n");
    fprintf(outfile, "\t---- ----- ------------------  ---------  -----------\n");
    for(h=0; h < moinfo.nirreps; h++) {
      for(i=0; i < params.rpi[h]; i++) {
        value = moinfo.uhf_evals[h][i];

        fprintf(outfile,   "\t%4d   %3s %18.14f  %9.5f  %10.2f\n", i, moinfo.labels[h],
            value, value*_hartree2ev, value*_hartree2wavenumbers);
      }
    }
    fprintf(outfile, "\n");

    fprintf(outfile, "\tUHF-CIS(D) Corrections:\n");
    fprintf(outfile, "\t-----------------------\n\n");

    fprintf(outfile, "\tRoot Irrep       Hartree          eV         cm-1\n");
    fprintf(outfile, "\t---- ----- ------------------  ---------  ----------\n");
    for(h=0; h < moinfo.nirreps; h++) {
      for(i=0; i < params.rpi[h]; i++) {
        value = moinfo.uhf_d[h][i];

        fprintf(outfile,   "\t%4d   %3s %18.14f  %9.5f  %10.2f\n", i, moinfo.labels[h],
            value, value*_hartree2ev, value*_hartree2wavenumbers);
      }
    }
    fprintf(outfile, "\n");

    fprintf(outfile, "\tUHF-CIS(D) Excitation Energies:\n");
    fprintf(outfile, "\t-------------------------------\n\n");

    fprintf(outfile, "\tRoot Irrep       Hartree          eV         cm-1\n");
    fprintf(outfile, "\t---- ----- ------------------  ---------  ----------\n");
    for(h=0; h < moinfo.nirreps; h++) {
      for(i=0; i < params.rpi[h]; i++) {
        value = moinfo.uhf_evals[h][i];
        value += moinfo.uhf_d[h][i];

        fprintf(outfile,   "\t%4d   %3s %18.14f  %9.5f  %10.2f\n", i, moinfo.labels[h],
            value, value*_hartree2ev, value*_hartree2wavenumbers);
      }
    }
    fprintf(outfile, "\n");

    free(evals);
    free(dcorr);
    free(symm);
  }

  dpd_close(0);
  if(params.ref == 0 || params.ref == 1) cachedone_rhf(cachelist);
  else if(params.ref == 2) cachedone_uhf(cachelist);
  free(cachefiles);

  cleanup();

  exit_io();
  exit(0);
}

void init_io(int argc, char *argv[])
{
  extern char *gprgid();
  char *progid;

  progid = (char *) malloc(strlen(gprgid())+2);
  sprintf(progid, ":%s",gprgid());

  psi_start(argc-1,argv+1,0);
  ip_cwk_add(":INPUT");
  ip_cwk_add(progid);
  free(progid);
  tstart(outfile);

  psio_init();

  psio_open(CC_INFO, 1);
  psio_open(CC_OEI, 1);
  psio_open(CC_CINTS, 1);
  psio_open(CC_DINTS, 1);
  psio_open(CC_EINTS, 1);
  psio_open(CC_FINTS, 1);
  psio_open(CC_DENOM, 1);
  psio_open(CC_MISC, 0);
  psio_open(CC_TMP0, 0);
  psio_open(CC_TMP1, 0);
}

void title(void)
{
  fprintf(outfile, "\t\t\t*************************\n");
  fprintf(outfile, "\t\t\t*                       *\n");
  fprintf(outfile, "\t\t\t*          CIS          *\n");
  fprintf(outfile, "\t\t\t*                       *\n");
  fprintf(outfile, "\t\t\t*************************\n");
}

void exit_io(void)
{
  psio_close(CC_INFO, 1);
  psio_close(CC_OEI, 1);
  psio_close(CC_CINTS, 1);
  psio_close(CC_DINTS, 1);
  psio_close(CC_EINTS, 1);
  psio_close(CC_FINTS, 1);
  psio_close(CC_DENOM, 1);
  psio_close(CC_MISC, 1);
  psio_close(CC_TMP0, 0);
  psio_close(CC_TMP1, 0);

  psio_done();
  tstop(outfile);
  psi_stop();
}

char *gprgid()
{
  char *prgid = "CIS";

  return(prgid);
}
