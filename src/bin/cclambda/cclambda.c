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
void init_io(void);
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
void update(void);
int converged(void);
int **cacheprep_rhf(int level, int *cachefiles);
int **cacheprep_uhf(int level, int *cachefiles);
void cachedone_rhf(int **cachelist);
void cachedone_uhf(int **cachelist);

int main(int argc, char *argv[])
{
  int done=0;
  int **cachelist, *cachefiles;

  moinfo.iter=0;
  
  init_io();
  title();
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

  params.maxiter = 1;
  for(moinfo.iter=1 ; moinfo.iter <= params.maxiter; moinfo.iter++) {
    sort_amps();
    G_build();
    L1_build();
    L2_build();
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
  overlap();
  dpd_close(0);
  cleanup(); 
  exit_io();
  exit(0);
}

void init_io(void)
{
  int i;
  char *gprgid();
  char *progid;

  progid = (char *) malloc(strlen(gprgid())+2);
  sprintf(progid, ":%s",gprgid());

  ffile(&infile,"input.dat",2);
  ffile(&outfile,"output.dat",1);
  tstart(outfile);
  ip_set_uppercase(1);
  ip_initialize(infile,outfile);
  ip_cwk_add(":DEFAULT");
  ip_cwk_add(progid);

  free(progid);

  psio_init();

  /* Open all dpd data files here */
  for(i=CC_MIN; i <= CC_MAX; i++) psio_open(i,1);
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
  ip_done();
  tstop(outfile);
  fclose(infile);
  fclose(outfile);
}

char *gprgid()
{
   char *prgid = "CCLAMBDA";

   return(prgid);
}
