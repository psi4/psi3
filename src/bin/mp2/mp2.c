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
#include <ccfiles.h>
#include "moinfo.h"
#include "params.h"
#include "globals.h"

void init_io(int argc, char *argv[]);
void title(void);
void get_moinfo(void);
void get_params(void);
void init_ioff(void);
int **cacheprep_rhf(int level, int *cachefiles);
void cachedone_rhf(int **cachelist);
struct dpd_file4_cache_entry *priority_list(void);
void energy(void);
void opdm(void);
void lag(void);
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

  cachelist = cacheprep_rhf(params.cachelev, cachefiles);

  priority = priority_list();
      
  dpd_init(0, mo.nirreps, params.memory, params.cachetype, cachefiles,
           cachelist, priority, 2, mo.actdoccpi, mo.actdoccsym,
	   mo.actvirtpi, mo.actvirtsym);
  
  energy();
  
  if (params.opdm) {
  opdm();
  }
  
  dpd_close(0);
  
  cleanup();

  exit_io();
  
  exit(0);
}

void init_io(int argc, char *argv[])
{
  int i;
  extern char *gprgid();
  char *progid;

  progid = (char *)malloc(strlen(gprgid())+2);
  sprintf(progid, ":%s",gprgid());

  psi_start(argc-1,argv+1,0);
  ip_cwk_add(progid);
  free(progid);
  tstart(outfile);
  
  psio_init();
  for(i=CC_MIN; i <= CC_LAMPS; i++) 
    psio_open(i,1);
  for(i=CC_HBAR; i <= CC_MAX; i++) 
    psio_open(i,0);
}

void title(void)
{
  fprintf(outfile, "\t\t\t*************************\n");
  fprintf(outfile, "\t\t\t*                       *\n");
  fprintf(outfile, "\t\t\t*          MP2          *\n");
  fprintf(outfile, "\t\t\t*                       *\n");
  fprintf(outfile, "\t\t\t*************************\n");
}

void init_ioff(void)
{
  int i;
  
  mo.ioff = init_int_array(MAXIOFF);
  mo.ioff[0] = 0;
  for(i=1; i < MAXIOFF; i++) {
    mo.ioff[i] = mo.ioff[i-1] + i;
  }
  
}

void cleanup(void)
{
  int i;
  
  free(params.wfn);
  free(params.ref);
  free(mo.fzdoccpi);
  free(mo.fzvirtpi);
  free(mo.actdoccpi);
  free(mo.actdoccsym);
  free(mo.actvirtpi);
  free(mo.actvirtsym);
  free(mo.doccpi);
  free(mo.virtpi);
  free(mo.mopi);
  for(i=0; i < mo.nirreps; i++)
    free(mo.irreplabels[i]);
  free(mo.irreplabels);
  free(mo.ioff);
  free(mo.scfevals);
}

void exit_io(void)
{
  int i;

  for(i=CC_MIN; i <= CC_MAX; i++) 
    psio_close(i,1);
  psio_done();
  tstop(outfile);
  psi_stop();
}

char *gprgid()
{
  char *prgid = "MP2";
  return(prgid);
}
