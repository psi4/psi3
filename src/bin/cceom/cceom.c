/*
**  CCEOM: Program to calculate the EOM CCSD right-hand eigenvector and  energy
*/
#include <stdio.h>
#include <stdlib.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "globals.h"

void init_io(void);
void get_moinfo(void);
void cleanup(void);
void exit_io(void);
void diag(void);
void hbar_clean();
void get_params(void);
void get_eom_params(void);
void form_dpd_dp(void);
int **cacheprep_uhf(int level, int *cachefiles);
int **cacheprep_rhf(int level, int *cachefiles);
void sort_amps(void);

/* local correlation functions */
void local_init(void);
void local_done(void);

int main(int argc, char *argv[])
{
  int i, h, done=0, *cachefiles, **cachelist;
  eom_params.dot_with_Lg = 0;
  for (i=1; i<argc; ++i) {
    if (!strcmp(argv[i],"--dot_with_Lg"))
      eom_params.dot_with_Lg = 1;
  }
  init_io();
  fprintf(outfile,"\n\t**********************************************************\n");
  fprintf(outfile,"\t*  CCEOM: An Equation of Motion Coupled Cluster Program  *\n");
  fprintf(outfile,"\t**********************************************************\n");

  get_moinfo();
  fflush(outfile);
  get_params();
  get_eom_params();
#ifdef TIME_CCEOM
  timer_init();
  timer_on("cceom");
#endif

  form_dpd_dp();
 
  cachefiles = init_int_array(PSIO_MAXUNIT);

  if (params.ref == 2) { // UHF 
    cachelist = cacheprep_uhf(params.cachelev, cachefiles);
    /* cachelist = init_int_matrix(32,32); */

    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles,
    cachelist, NULL, 4, moinfo.aoccpi, moinfo.aocc_sym, moinfo.avirtpi,
    moinfo.avir_sym, moinfo.boccpi, moinfo.bocc_sym, moinfo.bvirtpi, moinfo.bvir_sym);
  }
  else { // RHF or ROHF
    cachelist = cacheprep_rhf(params.cachelev, cachefiles);
    /* cachelist = init_int_matrix(12,12); */

    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles,
           cachelist, NULL, 2, moinfo.occpi, moinfo.occ_sym,
           moinfo.virtpi, moinfo.vir_sym);
  }

  if(params.local) local_init();

  diag();
  sort_amps();

  dpd_close(0);
  if(params.local) local_done();
  cleanup(); 
#ifdef TIME_CCEOM
  timer_off("cceom");
  timer_done();
#endif
  exit_io();
  exit(0);
}

void init_io(void)
{
  int i;
  char *gprgid();
  progid = (char *) malloc(strlen(gprgid())+2);
  sprintf(progid, ":%s",gprgid());
  ffile(&infile,"input.dat",2);
  ffile(&outfile,"output.dat",1);
  tstart(outfile);
  ip_set_uppercase(1);
  ip_initialize(infile,outfile);
  ip_cwk_clear();
  ip_cwk_add(":DEFAULT");
  ip_cwk_add(":INPUT");
  ip_cwk_add(progid);
  psio_init();
  for(i=CC_MIN; i <= CC_MISC; i++) psio_open(i,1);
  for(i=CC_TMP; i <= CC_MAX; i++) psio_open(i,0);
}

void exit_io(void)
{
  int i;
  for(i=CC_MIN; i <= CC_MISC; i++) psio_close(i,1);
  for(i=CC_TMP; i <= CC_MAX; i++) psio_close(i,1);
  psio_done();
  tstop(outfile);
  ip_done();
}

char *gprgid()
{
   char *prgid = "CCEOM";
   return(prgid);
}

void form_dpd_dp(void) {
  int h, h0, h1, cnt, nirreps;
  nirreps = moinfo.nirreps;

  dpd_dp = (int ***) malloc(nirreps * sizeof(int **));
  for(h=0; h < nirreps; h++) {
      dpd_dp[h] = init_int_matrix(nirreps,2);
      cnt=0;
      for(h0=0; h0 < nirreps; h0++) {
          for(h1=0; h1 < nirreps; h1++) {
              if((h0^h1)==h) {
                  dpd_dp[h][cnt][0] = h0;
                  dpd_dp[h][cnt++][1] = h1;
                }
            }
        }
    }
}

