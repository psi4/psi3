/*
**  CCRESPONSE: Program to compute CC linear response properties.
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
void init_io(int argc, char *argv[]);
void init_ioff(void);
void title(void);
void get_moinfo(void);
void get_params(void);
void cleanup(void);
void exit_io(void);
int **cacheprep_rhf(int level, int *cachefiles);
int **cacheprep_uhf(int level, int *cachefiles);
void cachedone_uhf(int **cachelist);
void cachedone_rhf(int **cachelist);
void transmu(void);
void sortmu(void);
void hbar_extra(void);
void sort_lamps(void);
void compute_X(char *cart, int irrep, double omega);
double LCX(char *cart, int irrep, double omega);
double HXY(char *cart, int irrep, double omega);
double LHX1Y1(char *cart, int irrep, double omega);
double LHX2Y2(char *cart, int irrep, double omega);
double LHX1Y2(char *cart, int irrep, double omega);
void mp2_density(void);

int main(int argc, char *argv[])
{
  int **cachelist, *cachefiles;
  double polar, polar_LCX, polar_HXY, polar_LHX1Y1, polar_LHX2Y2, polar_LHX1Y2;

  init_io(argc, argv);
  init_ioff();
  title();
  get_moinfo();
  get_params();

  cachefiles = init_int_array(PSIO_MAXUNIT);

  if(params.ref == 2) { /*** UHF references ***/
    cachelist = cacheprep_uhf(params.cachelev, cachefiles);

    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, 
	     NULL, 4, moinfo.aoccpi, moinfo.aocc_sym, moinfo.avirtpi, moinfo.avir_sym,
	     moinfo.boccpi, moinfo.bocc_sym, moinfo.bvirtpi, moinfo.bvir_sym);
  } 
  else { /*** RHF/ROHF references ***/
    cachelist = cacheprep_rhf(params.cachelev, cachefiles);

    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL,
	     2, moinfo.occpi, moinfo.occ_sym, moinfo.virtpi, moinfo.vir_sym);
  }

  transmu();
  sortmu();
  mubar();
  hbar_extra();
  sort_lamps();

  polar = 0.0;

  compute_X("X", moinfo.irrep_x, params.omega);
  if(params.omega != 0.0) 
    compute_X("X", moinfo.irrep_x, -params.omega);

  polar_LCX = LCX("X", moinfo.irrep_x, params.omega);
  polar_LCX += LCX("X", moinfo.irrep_x, -params.omega);
  fprintf(outfile, "polar_LCX = %20.12f\n", polar_LCX);

  polar += polar_LCX;

  polar_HXY = HXY("X", moinfo.irrep_x, params.omega);
  fprintf(outfile, "polar_HXY = %20.12f\n", polar_HXY);

  polar += polar_HXY;

  polar_LHX1Y1 = LHX1Y1("X", moinfo.irrep_x, params.omega);
  fprintf(outfile, "polar_LHX1Y1 = %20.12f\n", polar_LHX1Y1);

  polar += polar_LHX1Y1;

  polar_LHX2Y2 = LHX2Y2("X", moinfo.irrep_x, params.omega);
  fprintf(outfile, "polar_LHX2Y2 = %20.12f\n", polar_LHX2Y2);

  polar += polar_LHX2Y2;

  polar_LHX1Y2 = LHX1Y2("X", moinfo.irrep_x, params.omega);
  polar_LHX1Y2 += LHX1Y2("X", moinfo.irrep_x, -params.omega);
  fprintf(outfile, "polar_LHX1Y2 = %20.12f\n", polar_LHX1Y2);

  polar += polar_LHX1Y2;
  fprintf(outfile, "%1s%1s polar = %20.12f\n", "X", "X", polar);

  polar = 0.0;

  compute_X("Y", moinfo.irrep_y, params.omega);
  if(params.omega != 0.0) 
    compute_X("Y", moinfo.irrep_y, -params.omega);

  polar_LCX = LCX("Y", moinfo.irrep_y, params.omega);
  polar_LCX += LCX("Y", moinfo.irrep_y, -params.omega);
  fprintf(outfile, "polar_LCX = %20.12f\n", polar_LCX);

  polar += polar_LCX;

  polar_HXY = HXY("Y", moinfo.irrep_y, params.omega);
  fprintf(outfile, "polar_HXY = %20.12f\n", polar_HXY);

  polar += polar_HXY;

  polar_LHX1Y1 = LHX1Y1("Y", moinfo.irrep_y, params.omega);
  fprintf(outfile, "polar_LHX1Y1 = %20.12f\n", polar_LHX1Y1);

  polar += polar_LHX1Y1;

  polar_LHX2Y2 = LHX2Y2("Y", moinfo.irrep_y, params.omega);
  fprintf(outfile, "polar_LHX2Y2 = %20.12f\n", polar_LHX2Y2);

  polar += polar_LHX2Y2;

  polar_LHX1Y2 = LHX1Y2("Y", moinfo.irrep_y, params.omega);
  polar_LHX1Y2 += LHX1Y2("Y", moinfo.irrep_y, -params.omega);
  fprintf(outfile, "polar_LHX1Y2 = %20.12f\n", polar_LHX1Y2);

  polar += polar_LHX1Y2;
  fprintf(outfile, "%1s%1s polar = %20.12f\n", "Y", "Y", polar);

  polar = 0.0;

  compute_X("Z", moinfo.irrep_z, params.omega);
  if(params.omega != 0.0) 
    compute_X("Z", moinfo.irrep_z, -params.omega);

  polar_LCX = LCX("Z", moinfo.irrep_z, params.omega);
  polar_LCX += LCX("Z", moinfo.irrep_z, -params.omega);
  fprintf(outfile, "polar_LCX = %20.12f\n", polar_LCX);

  polar += polar_LCX;

  polar_HXY = HXY("Z", moinfo.irrep_z, params.omega);
  fprintf(outfile, "polar_HXY = %20.12f\n", polar_HXY);

  polar += polar_HXY;

  polar_LHX1Y1 = LHX1Y1("Z", moinfo.irrep_z, params.omega);
  fprintf(outfile, "polar_LHX1Y1 = %20.12f\n", polar_LHX1Y1);

  polar += polar_LHX1Y1;

  polar_LHX2Y2 = LHX2Y2("Z", moinfo.irrep_z, params.omega);
  fprintf(outfile, "polar_LHX2Y2 = %20.12f\n", polar_LHX2Y2);

  polar += polar_LHX2Y2;

  polar_LHX1Y2 = LHX1Y2("Z", moinfo.irrep_z, params.omega);
  polar_LHX1Y2 += LHX1Y2("Z", moinfo.irrep_z, -params.omega);
  fprintf(outfile, "polar_LHX1Y2 = %20.12f\n", polar_LHX1Y2);

  polar += polar_LHX1Y2;
  fprintf(outfile, "%1s%1s polar = %20.12f\n", "Z", "Z", polar);

  dpd_close(0);

  if(params.ref == 2) cachedone_uhf(cachelist);
  else cachedone_rhf(cachelist);
  free(cachefiles);

  cleanup();

  exit_io();
  exit(PSI_RETURN_SUCCESS);
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

  for(i=CC_MIN; i <= CC_MAX; i++) psio_open(i, 1);
}

void title(void)
{
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t*       CCRESPONSE       *\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t**************************\n");
}

void exit_io(void)
{
  int i;

  for(i=CC_MIN; i <= CC_MAX; i++) psio_close(i, 1);

  psio_done();
  tstop(outfile);
  psi_stop();
}

char *gprgid()
{
   char *prgid = "CCRESPONSE";

   return(prgid);
}

void init_ioff(void)
{
  int i;
  ioff = init_int_array(IOFF_MAX);
  ioff[0] = 0;
  for(i=1; i < IOFF_MAX; i++) ioff[i] = ioff[i-1] + i;
}
