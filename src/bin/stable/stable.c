/*
**  STABLE: Program to calculate the eigenvalues and eigenvectors
**  of the molecular orbital Hessian.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libchkpt/chkpt.h>
#include <psifiles.h>
#include "globals.h"

/* Max length of ioff array */
#define IOFF_MAX 32641

/* Function prototypes */
void init_io(void);
void init_ioff(void);
void title(void);
void get_moinfo(void);
void get_params(void);
void exit_io(void);
void sort_oei(void);
void sort_tei(void);
void fock(void);
void build_A(void);
void diag_A(void);
int **cacheprep(int level, int *cachefiles);

int main(int argc, char *argv[])
{
  int **cachelist, *cachefiles;

  init_io(argc, argv);
  init_ioff();
  title();
  get_moinfo();
  get_params();

  cachefiles = init_int_array(PSIO_MAXUNIT);
  cachelist = cacheprep(params.cachelev, cachefiles);
  dpd_init(0, moinfo.nirreps, params.memory, cachefiles, cachelist,
           2, moinfo.occpi, moinfo.occ_sym, moinfo.virtpi, moinfo.vir_sym);

  sort_oei();
  sort_tei();
  fock();
  build_A();
  diag_A();

  dpd_close(0);
  exit_io();
  exit(0);
}

void init_io(int argc, char *argv[])
{
  int i;
  extern char *gprgid();
  char *progid;

  progid = (char *) malloc(strlen(gprgid())+2);
  sprintf(progid, ":%s",gprgid());

  psi_start(argc-1,argv+1,0); /* this assumes no cmdline args except filenames */
  ip_cwk_add(progid);
  free(progid);
  tstart(outfile);
  psio_init();

  psio_open(PSIF_MO_HESS,0);
}

void title(void)
{
  fprintf(outfile, "\n");
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t*         STABLE         *\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\n");
}

void exit_io(void)
{
  int i;
 
  psio_close(PSIF_MO_HESS,1);

  psio_done();
  tstop(outfile);
  psi_stop();
}

char *gprgid()
{
   char *prgid = "STABLE";

   return(prgid);
}

void init_ioff(void)
{
  int i;
  ioff = init_int_array(IOFF_MAX);
  ioff[0] = 0;
  for(i=1; i < IOFF_MAX; i++) ioff[i] = ioff[i-1] + i;
}
