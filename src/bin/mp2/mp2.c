#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#include <libiwl/iwl.h>
#include <physconst.h>
#include <psifiles.h>
#include "globals.h"

void init_io(int argc, char *argv[]);
void title(void);
void get_moinfo(void);
void get_params(void);
void init_ioff(void);
int **cacheprep_rhf(int level, int *cachefiles);
int **cacheprep_uhf(int level, int *cachefiles);
void cachedone_rhf(int **cachelist);
struct dpd_file4_cache_entry *priority_list(void);
double energy(void);
void sort_amps(void);
void opdm(void);
void twopdm(void);
void lag(void);
void build_X(void);
void build_A(void);
void Zvector(void);
void relax_I(void);
void relax_opdm(void);
void sort_opdm(void);
void sort_I(void);
void fold(void);
void deanti(void);
void write_data(void);
void check_energy(int);
void sort_twopdm(void);
void cleanup(void);
void exit_io(void);

int main(int argc, char *argv[])
{
  int *cachefiles;
  int **cachelist;
 
  struct dpd_file4_cache_entry *priority;
  
  init_io(argc,argv);
  title();

  get_moinfo();
  get_params();
  init_ioff();
  
  cachefiles = init_int_array(PSIO_MAXUNIT);

  if(params.ref == 2) { /** UHF **/
    cachelist = cacheprep_uhf(params.cachelev,cachefiles);
    dpd_init(0,mo.nirreps,params.memory,params.cachetype,cachefiles,cachelist,
             NULL,4,mo.aoccpi,mo.aocc_sym,mo.avirpi,mo.avir_sym,mo.boccpi,
             mo.bocc_sym,mo.bvirpi,mo.bvir_sym);
  }
  else { /** RHF or ROHF **/
    cachelist = cacheprep_rhf(params.cachelev,cachefiles);
    priority = priority_list();
    dpd_init(0,mo.nirreps,params.memory,params.cachetype,cachefiles,
             cachelist,priority,2,mo.occpi,mo.occ_sym,mo.virpi,mo.vir_sym);
  }
  
  amps();

  mo.Emp2 = energy();
  
  fprintf(outfile,"\n");
  fprintf(outfile,"\tMP2 correlation energy  = %20.15f\n",mo.Emp2);
  fprintf(outfile,"\tMP2 total energy        = %20.15f\n",mo.Escf+mo.Emp2);
  fflush(outfile);

  chkpt_init(PSIO_OPEN_OLD);
  chkpt_wt_etot(mo.Escf+mo.Emp2);
  chkpt_close();
  
  if(params.opdm) {
    opdm();
    if(params.relax_opdm) {
      lag();
      build_A();
      Zvector();
    }
    sort_opdm();
    //dipole();
  }

  if(params.gradient) {
    sort_amps();
    opdm();
    twopdm();
    lag();
    build_X();
    build_A();
    Zvector();
    relax_I();
    relax_opdm();
    sort_I(); 
    sort_opdm();
    fold();
    deanti();
    write_data();
  }

  dpd_close(0);
  
  cleanup();

  exit_io();
  
  exit(0);
}

void init_io(int argc, char *argv[])
{
  int i=0;
  int num_extra_args=0;
  char **extra_args;
  extern char *gprgid();
  char *progid;

  extra_args = (char **)malloc(argc*sizeof(char *));
  progid = (char *)malloc(strlen(gprgid())+2);
  sprintf(progid, ":%s",gprgid());
 
  params.opdm = 0;

  for(i=1; i<argc; i++) {
    if(strcmp(argv[i], "--opdm") == 0) {
      params.opdm = 1;
    }
    else {
      extra_args[num_extra_args++] = argv[i];
    }
  }
  
  psi_start(num_extra_args,extra_args,0);
  ip_cwk_add(progid);
  free(progid);
  tstart(outfile);
  
  psio_init();
  for(i=CC_MIN; i <= CC_MAX; i++) 
    psio_open(i,1);

  free(extra_args);
}

void title(void)
{
  fprintf(outfile, "\t\t\t*************************\n");
  fprintf(outfile, "\t\t\t*                       *\n");
  fprintf(outfile, "\t\t\t*          MP2          *\n");
  fprintf(outfile, "\t\t\t*                       *\n");
  fprintf(outfile, "\t\t\t*************************\n");
  fflush(outfile);
}

void init_ioff(void)
{
  int i;
  
  ioff = init_int_array(MAXIOFF);
  ioff[0] = 0;
  for(i=1; i < MAXIOFF; i++) {
    ioff[i] = ioff[i-1] + i;
  }
  
}

void cleanup(void)
{
  int i;
  
  free(params.wfn);
  
  free(mo.doccpi);
  free(mo.soccpi);
  free(mo.mopi);
  free(mo.fzdoccpi);
  free(mo.fzvirtpi);
  for(i=0; i < mo.nirreps; i++)
    free(mo.irreplabels[i]);
  free(mo.irreplabels);
  free(ioff);
  
  if(params.ref == 2) {
    free(mo.aoccpi);
    free(mo.boccpi);
    free(mo.avirpi);
    free(mo.bvirpi);
    free(mo.aocc_sym);
    free(mo.bocc_sym);
    free(mo.avir_sym);
    free(mo.bvir_sym);
    free(mo.aocc_off);
    free(mo.bocc_off);
    free(mo.avir_off);
    free(mo.bvir_off);
    free(mo.qt_aocc);
    free(mo.qt_bocc);
    free(mo.qt_avir);
    free(mo.qt_bvir);
  }
  else {
    free(mo.occpi);
    free(mo.virpi);
    free(mo.occ_sym);
    free(mo.vir_sym);
    free(mo.occ_off);
    free(mo.vir_off);
    free(mo.qt_occ);
    free(mo.qt_vir);
  }
}

void exit_io(void)
{
  int i;

  for(i=CC_MIN; i <= CC_MAX; i++) 
    psio_close(i,0);
  psio_done();
  tstop(outfile);
  psi_stop();
}

char *gprgid()
{
  char *prgid = "MP2";
  return(prgid);
}
