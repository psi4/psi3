/*
**  CCLAMBDA: Program to calculate the coupled-cluster lambda vector.
*/

#include <stdio.h>
#include <stdlib.h>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include "globals.h"

/* Function prototypes */
void init_io(int argc, char *argv[]);
void title(void);
void get_moinfo(void);
void get_params(void);
void cleanup(void);
void init_amps(void);
double pseudoenergy(void);
void exit_io(void);
void G_build();
void L1_build(void);
void L2_build(void);
void sort_amps(void);
void Lsave(void);
void Lnorm(void);
void Lmag(void);
void update(void);
int converged(void);
int **cacheprep_rhf(int level, int *cachefiles);
int **cacheprep_uhf(int level, int *cachefiles);
void cachedone_rhf(int **cachelist);
void cachedone_uhf(int **cachelist);

int main(int argc, char *argv[])
{
  int done=0, i;
  int **cachelist, *cachefiles;

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

  init_amps();
  fprintf(outfile, "\t          Solving Lambda Equations\n");
  fprintf(outfile, "\t          ------------------------\n");
  fprintf(outfile, "\tIter         PseudoEnergy             RMS  \n");
  fprintf(outfile, "\t----     ---------------------     --------\n");
  moinfo.lcc = pseudoenergy();
  update();

  denom();

  for(moinfo.iter=1 ; moinfo.iter <= params.maxiter; moinfo.iter++) {
    sort_amps();
    G_build();
    L1_build();
    L2_build();
    if (!params.ground) Lmag();
    if(converged()) {
      done = 1;  /* Boolean for convergence */
      Lsave();
      moinfo.lcc = pseudoenergy();
      sort_amps();  /* For later calculations */
      update();
      fprintf(outfile, "\n\tIterations converged.\n");
      fflush(outfile);
      break;
    }
/*
    if (!params.ground) Lnorm();
*/
    diis(moinfo.iter);
    Lsave();
    moinfo.lcc = pseudoenergy();
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
    exit(1);
  }
  if (params.ground) overlap();
  if (!params.ground) Lnorm();
  dpd_close(0);

  psio_write_entry(CC_INFO,"EOM L0",(char *) &params.L0, sizeof(double));
  psio_write_entry(CC_INFO,"EOM L Irrep",(char *) &L_irr, sizeof(int));

  fprintf(outfile,"Writing EOM L0 %15.10lf\n", params.L0);

  if(params.ref == 2) cachedone_uhf(cachelist);
  else cachedone_rhf(cachelist);
  free(cachefiles);

  cleanup(); 
  exit_io();
  exit(0);
}

/* parse command line arguments */
void init_io(int argc, char *argv[])
{
  int i, num_unparsed;
  extern char *gprgid();
  char *progid, *argv_unparsed[100];

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
  ip_cwk_add(progid);
  free(progid);
  tstart(outfile);
  psio_init();

  for(i=CC_MIN; i <= CC_MAX; i++) psio_open(i,1);

  params.cceom_energy = 0.0;
  if (params.ground) {
    L_irr = 0;
    params.L0 = 1.0;
    params.cceom_energy = 0.0;
  }
  else {
    /* assume symmetry of L is that of R */
    psio_read_entry(CC_INFO,"EOM R Irrep", (char *) &L_irr, sizeof(int));
    params.L0 = 0.0;
    psio_read_entry(CC_INFO,"EOM R0", (char *) &(params.R0),sizeof(double));
    psio_read_entry(CC_INFO,"CCEOM Energy",
      (char *) &(params.cceom_energy),sizeof(double));
  }
}

void title(void)
{
  fprintf(outfile, "\n");
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t*        CCLAMBDA        *\n");
  fprintf(outfile, "\t\t\t*                        *\n");
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
