/*
**  CCLAMBDA: Program to calculate the coupled-cluster lambda vector.
*/

#include <stdio.h>
#include <stdlib.h>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <psifiles.h>
#include "globals.h"

/* Function prototypes */
void init_io(int argc, char *argv[]);
void title(void);
void get_moinfo(void);
void get_params(void);
void cleanup(void);
void init_amps(int L_irr, int state__index);
double pseudoenergy(int L_irr);
void exit_io(void);
void G_build(int L_irr);
void L1_build(int L_irr, int root_L_irr);
void L2_build(int L_irr, int root_L_irr);
void sort_amps(int L_irr);
void Lsave(int L_irr);
void Lnorm(int L_irr, int root_L_irr);
void Lmag(void);
void update(void);
int converged(int L_irr);
int **cacheprep_rhf(int level, int *cachefiles);
int **cacheprep_uhf(int level, int *cachefiles);
void cachedone_rhf(int **cachelist);
void cachedone_uhf(int **cachelist);
void denom(int L_irr, int e_index);
void overlap(int L_irr);
void Lsave_index(int L_irr, int root_L_irr);
void L_to_LAMPS(int L_irr);
extern void check_ortho(void);
void L_zero(int L_irr);

int main(int argc, char *argv[])
{
  int done=0, i, L_irr, root_L_irr;
  int **cachelist, *cachefiles;
  dpdfile2 L1;

  init_io(argc, argv); /* parses command-line arguments */
  title();
  moinfo.iter=0;
  get_moinfo();
  get_params();

  cachefiles = init_int_array(PSIO_MAXUNIT);

  if(params.ref == 0 || params.ref == 1) { /** RHF or ROHF **/

    cachelist = cacheprep_rhf(params.cachelev, cachefiles);

    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL,
	     2, moinfo.occpi, moinfo.occ_sym, moinfo.virtpi, moinfo.vir_sym);

    if(params.aobasis) { /* Set up new DPD for AO-basis algorithm */
      dpd_init(1, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL, 
	       2, moinfo.occpi, moinfo.orbsym, moinfo.orbspi, moinfo.orbsym);
      dpd_set_default(0);
    }

  }
  else if(params.ref == 2) { /** UHF **/

    cachelist = cacheprep_uhf(params.cachelev, cachefiles);

    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, 
	     cachelist, NULL, 4, moinfo.aoccpi, moinfo.aocc_sym, moinfo.avirtpi,
	     moinfo.avir_sym, moinfo.boccpi, moinfo.bocc_sym, moinfo.bvirtpi, moinfo.bvir_sym);

    if(params.aobasis) { /* Set up new DPD's for AO-basis algorithm */
      dpd_init(1, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL, 
               4, moinfo.aoccpi, moinfo.aocc_sym, moinfo.orbspi, moinfo.orbsym, 
	       moinfo.boccpi, moinfo.bocc_sym, moinfo.orbspi, moinfo.orbsym);
      dpd_set_default(0);
    }
  }


  if(params.local) local_init();

  for (L_irr=0; L_irr<moinfo.nirreps; ++L_irr) {
    psio_close(CC_TMP,0);
    psio_close(CC_TMP1,0);
    psio_close(CC_TMP2,0);
    psio_close(CC_LAMBDA,0);
    psio_open(CC_TMP,0);
    psio_open(CC_TMP1,0);
    psio_open(CC_TMP2,0);
    psio_open(CC_LAMBDA,0);
    for (root_L_irr=0; root_L_irr < params.Ls_per_irrep[L_irr]; ++root_L_irr) {
      fprintf(outfile,"\tSymmetry of excited state: %s\n", moinfo.labels[moinfo.sym^L_irr]);
      fprintf(outfile,"\tSymmetry of right eigenvector: %s\n",moinfo.labels[L_irr]);

      init_amps(L_irr, root_L_irr);

      fprintf(outfile, "\n\t          Solving Lambda Equations\n");
      fprintf(outfile, "\t          ------------------------\n");
      fprintf(outfile, "\tIter     PseudoEnergy or Norm         RMS  \n");
      fprintf(outfile, "\t----     ---------------------     --------\n");

      denom(L_irr, root_L_irr); /* second argument determines E used */

      moinfo.lcc = pseudoenergy(L_irr);
      update();

      for(moinfo.iter=1 ; moinfo.iter <= params.maxiter; moinfo.iter++) {
        sort_amps(L_irr);
        G_build(L_irr);

        /* must zero New L before adding RHS */
        L_zero(L_irr);
        L1_build(L_irr,root_L_irr);
        L2_build(L_irr,root_L_irr);
  
        if(converged(L_irr)) {
          done = 1;  /* Boolean for convergence */
          Lsave(L_irr); /* copy "New L" to "L" */
          moinfo.lcc = pseudoenergy(L_irr);
          update();
          if (!params.ground) {
            Lnorm(L_irr, root_L_irr); /* normalize against R */
            Lsave_index(L_irr, root_L_irr); /* put Ls in unique location */
          }
          else {
            L_to_LAMPS(L_irr);
          }
         /* sort_amps(); to be done by later functions */
          fprintf(outfile, "\n\tIterations converged.\n");
          fflush(outfile);
          break;
        }
  
        diis(moinfo.iter, L_irr);
        Lsave(L_irr);
        moinfo.lcc = pseudoenergy(L_irr);
        update();
      }
      fprintf(outfile, "\n");
      if(!done) {
        fprintf(outfile, "\t ** Lambda not converged to %2.1e ** \n",
            params.convergence);
        fflush(outfile);
        dpd_close(0);
        cleanup();
        exit_io();
        exit(PSI_RETURN_FAILURE);
      }
      if (params.ground) {
        overlap(L_irr);
        overlap_LAMPS(L_irr);
      }
    }
  }

  if(params.local) local_done();

  if (!params.ground) check_ortho();

  dpd_close(0);

  if(params.ref == 2) cachedone_uhf(cachelist);
  else cachedone_rhf(cachelist);
  free(cachefiles);

  cleanup(); 
  exit_io();
  exit(PSI_RETURN_SUCCESS);
}

/* parse command line arguments */
void init_io(int argc, char *argv[])
{
  int i, num_unparsed;
  extern char *gprgid();
  char *lbl, *progid, *argv_unparsed[100];

  progid = (char *) malloc(strlen(gprgid())+2);
  sprintf(progid, ":%s",gprgid());

  params.ground = 1;
  for (i=1, num_unparsed=0; i<argc; ++i) {
    if (!strcmp(argv[i],"--excited")) {
      params.ground = 0;
    }
    else {
      argv_unparsed[num_unparsed++] = argv[i];
    }
  }

  psi_start(num_unparsed, argv_unparsed, 0);
  ip_cwk_add(":INPUT");
  ip_cwk_add(progid);
  free(progid);

  ip_cwk_add(":CCEOM");

  tstart(outfile);
  psio_init();

  for(i=CC_MIN; i <= CC_MAX; i++) psio_open(i,1);

  if (params.ground)
    params.L0 = 1.0;
  else
    params.L0 = 0.0;

}

void title(void)
{
  fprintf(outfile, "\n");
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\t\t\t*        CCLAMBDA        *\n");
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\n");
}

void exit_io(void)
{
  int i;
 
  /* Close all dpd data files here */
  for(i=CC_MIN; i <= CC_MAX; i++) psio_close(i,1);

  psio_done();
  tstop(outfile);
  psi_stop();
}

char *gprgid()
{
   char *prgid = "CCLAMBDA";

   return(prgid);
}

/* put copies of L for excited states in LAMPS with irrep and index label */
void Lsave_index(int L_irr, int root_L_irr) {
  dpdfile2 L1;
  dpdbuf4 L2, LIjAb, LIjbA;
  char L1A_lbl[32], L1B_lbl[32], L2AA_lbl[32], L2BB_lbl[32], L2AB_lbl[32];
  char lbl[32];

  sprintf(L1A_lbl, "LIA %d %d", L_irr, root_L_irr);
  sprintf(L1B_lbl, "Lia %d %d", L_irr, root_L_irr);
  sprintf(L2AA_lbl, "LIJAB %d %d", L_irr, root_L_irr);
  sprintf(L2BB_lbl, "Lijab %d %d", L_irr, root_L_irr);
  sprintf(L2AB_lbl, "LIjAb %d %d", L_irr, root_L_irr);

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/
    dpd_file2_init(&L1, CC_OEI, L_irr, 0, 1, "LIA");
    dpd_file2_copy(&L1, CC_OEI, L1A_lbl);
    dpd_file2_close(&L1);
    dpd_file2_init(&L1, CC_OEI, L_irr, 0, 1, "Lia");
    dpd_file2_copy(&L1, CC_OEI, L1B_lbl);
    dpd_file2_close(&L1);
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_copy(&L2, CC_LAMPS, L2AA_lbl);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
    dpd_buf4_copy(&L2, CC_LAMPS, L2BB_lbl);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_copy(&L2, CC_LAMPS, L2AB_lbl);
    dpd_buf4_close(&L2);
  }
  else if(params.ref == 2) { /** UHF **/
    dpd_file2_init(&L1, CC_OEI, L_irr, 0, 1, "LIA");
    dpd_file2_copy(&L1, CC_OEI, L1A_lbl);
    dpd_file2_close(&L1);
    dpd_file2_init(&L1, CC_OEI, L_irr, 2, 3, "Lia");
    dpd_file2_copy(&L1, CC_OEI, L1B_lbl);
    dpd_file2_close(&L1);
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_copy(&L2, CC_LAMPS, L2AA_lbl);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "Lijab");
    dpd_buf4_copy(&L2, CC_LAMPS, L2BB_lbl);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_copy(&L2, CC_LAMPS, L2AB_lbl);
    dpd_buf4_close(&L2);
  }

  if (params.ref == 0) { /** RHF for those codes that can use them **/
    dpd_buf4_init(&LIjAb, CC_LAMPS, L_irr, 0, 5, 0, 5, 0, L2AB_lbl);
    dpd_buf4_sort(&LIjAb, CC_TMP, pqsr, 0, 5, "LIjbA");
    sprintf(lbl, "2LIjAb - LIjbA %d %d", L_irr, root_L_irr);
    dpd_buf4_copy(&LIjAb, CC_LAMPS, lbl);
    dpd_buf4_close(&LIjAb);

    sprintf(lbl, "2LIjAb - LIjbA %d %d", L_irr, root_L_irr);
    dpd_buf4_init(&LIjAb, CC_LAMPS, L_irr, 0, 5, 0, 5, 0, lbl);
    dpd_buf4_scm(&LIjAb, 2.0);
    dpd_buf4_init(&LIjbA, CC_TMP, L_irr, 0, 5, 0, 5, 0, "LIjbA");
    dpd_buf4_axpy(&LIjbA, &LIjAb, -1.0);
    dpd_buf4_close(&LIjbA);
    dpd_buf4_close(&LIjAb);
  }

/* also put copy of L1 in LAMPS - RAK thinks this is better organization */
  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/
    dpd_file2_init(&L1, CC_OEI, L_irr, 0, 1, "LIA");
    dpd_file2_copy(&L1, CC_LAMPS, L1A_lbl); 
    dpd_file2_close(&L1);
    dpd_file2_init(&L1, CC_OEI, L_irr, 0, 1, "Lia");
    dpd_file2_copy(&L1, CC_LAMPS, L1B_lbl);
    dpd_file2_close(&L1);
  }
  else if(params.ref == 2) { /** UHF **/
    dpd_file2_init(&L1, CC_OEI, L_irr, 0, 1, "LIA");
    dpd_file2_copy(&L1, CC_LAMPS, L1A_lbl);
    dpd_file2_close(&L1);
    dpd_file2_init(&L1, CC_OEI, L_irr, 2, 3, "Lia");
    dpd_file2_copy(&L1, CC_LAMPS, L1B_lbl);
    dpd_file2_close(&L1);
  }
  return;
}

/* just used for ground state */
void L_to_LAMPS(int L_irr) {
  dpdfile2 L1;
  dpdbuf4 L2, LIjAb, LIjbA;

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/
    dpd_file2_init(&L1, CC_OEI, L_irr, 0, 1, "LIA");
    dpd_file2_copy(&L1, CC_LAMPS, "LIA");
    dpd_file2_close(&L1);
    dpd_file2_init(&L1, CC_OEI, L_irr, 0, 1, "Lia");
    dpd_file2_copy(&L1, CC_LAMPS, "Lia");
    dpd_file2_close(&L1);
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_copy(&L2, CC_LAMPS, "LIJAB");
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
    dpd_buf4_copy(&L2, CC_LAMPS, "Lijab");
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_copy(&L2, CC_LAMPS, "LIjAb");
    dpd_buf4_close(&L2);
  }
  else if(params.ref == 2) { /** UHF **/
    dpd_file2_init(&L1, CC_OEI, L_irr, 0, 1, "LIA");
    dpd_file2_copy(&L1, CC_OEI, "LIA");
    dpd_file2_close(&L1);
    dpd_file2_init(&L1, CC_OEI, L_irr, 2, 3, "Lia");
    dpd_file2_copy(&L1, CC_OEI, "Lia");
    dpd_file2_close(&L1);
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_copy(&L2, CC_LAMPS, "LIJAB");
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "Lijab");
    dpd_buf4_copy(&L2, CC_LAMPS, "Lijab");
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_copy(&L2, CC_LAMPS, "LIjAb");
    dpd_buf4_close(&L2);
  }
  if (params.ref == 0) { /** RHF for those codes that can use them **/
    dpd_buf4_init(&L2, CC_LAMPS, 0, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_scmcopy(&L2, CC_LAMPS, "2 LIjAb - LIjBa", 2);
    dpd_buf4_sort_axpy(&L2, CC_LAMPS, pqsr, 0, 5, "2 LIjAb - LIjBa", -1);
    dpd_buf4_close(&L2);

    /*
    dpd_buf4_init(&LIjAb, CC_LAMPS, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_sort(&LIjAb, CC_TMP, pqsr, 0, 5, "LIjbA");
    dpd_buf4_copy(&LIjAb, CC_LAMPS, "2LIjAb - LIjbA");
    dpd_buf4_close(&LIjAb);
    dpd_buf4_init(&LIjAb, CC_LAMPS, L_irr, 0, 5, 0, 5, 0, "2LIjAb - LIjbA");
    dpd_buf4_scm(&LIjAb, 2.0);
    dpd_buf4_init(&LIjbA, CC_TMP, L_irr, 0, 5, 0, 5, 0, "LIjbA");
    dpd_buf4_axpy(&LIjbA, &LIjAb, -1.0);
    dpd_buf4_close(&LIjbA);
    dpd_buf4_close(&LIjAb);
    */
  }
  return;
}

void L_zero(int L_irr) {
  dpdfile2 LIA, Lia;
  dpdbuf4 LIJAB, Lijab, LIjAb;

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/
    dpd_file2_init(&LIA, CC_OEI, L_irr, 0, 1, "New LIA");
    dpd_file2_init(&Lia, CC_OEI, L_irr, 0, 1, "New Lia");
  }
  else if(params.ref == 2) { /** UHF **/
    dpd_file2_init(&LIA, CC_OEI, L_irr, 0, 1, "New LIA");
    dpd_file2_init(&Lia, CC_OEI, L_irr, 2, 3, "New Lia");
  }
  dpd_file2_scm(&LIA, 0.0);
  dpd_file2_scm(&Lia, 0.0);
  dpd_file2_close(&LIA);
  dpd_file2_close(&Lia);

  if (params.ref == 0 || params.ref == 1 ) { /** RHF/ROHF **/
    dpd_buf4_init(&LIJAB, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    dpd_buf4_scm(&LIJAB, 0.0);
    dpd_buf4_close(&LIJAB);
    dpd_buf4_init(&Lijab, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab");
    dpd_buf4_scm(&Lijab, 0.0);
    dpd_buf4_close(&Lijab);
    dpd_buf4_init(&LIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    dpd_buf4_scm(&LIjAb, 0.0);
    dpd_buf4_close(&LIjAb);
  }
  else { /** UHF **/
    dpd_buf4_init(&LIJAB, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    dpd_buf4_scm(&LIJAB, 0.0);
    dpd_buf4_close(&LIJAB);
    dpd_buf4_init(&Lijab, CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "New Lijab");
    dpd_buf4_scm(&Lijab, 0.0);
    dpd_buf4_close(&Lijab);
    dpd_buf4_init(&LIjAb, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "New LIjAb");
    dpd_buf4_scm(&LIjAb, 0.0);
    dpd_buf4_close(&LIjAb);
  }
}
