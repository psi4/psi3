/*
**  CCSORT: Program to reorganize integrals for CC and MBPT calculations.
*/

#include <stdio.h>
#include <stdlib.h>
#include <ip_libv1.h>
#include <libciomr.h>
#include <dpd.h>
#include <psio.h>
#include <qt.h>
#include "MOInfo.h"
#include "Params.h"
#include "globals.h"

/* Max length of ioff array */
#define IOFF_MAX 32641

/* Function prototypes */
void init_io(void);
void init_ioff(void);
void title(void);
void get_params(void);
void get_moinfo(void);
void sort_oei(void);
void sort_tei(void);
void b_sort(void);
void c_sort(void);
void d_sort(void);
void e_sort(void);
void f_sort(void);
void scf_check(void);
void fock(void);
void denom(void);
void exit_io(void);
void cleanup(void);
int **cacheprep(int level, int *cachefiles);

int main(int argc, char *argv[])
{
  int **cachelist, *cachefiles;

  init_io();
  init_ioff();
  title();
  timer_init();
  get_moinfo();
  get_params();
  cachefiles = init_int_array(PSIO_MAXUNIT);
  cachelist = cacheprep(params.cachelev, cachefiles);
  dpd_init(0, moinfo.nirreps, params.memory, cachefiles, cachelist, 
           2, moinfo.occpi, moinfo.occ_sym, moinfo.virtpi, moinfo.vir_sym);
  sort_oei();
  sort_tei();
  b_sort();
  c_sort();
  d_sort();
  e_sort();
  f_sort();
  scf_check();
  fock();
  denom();
  dpd_close(0);
  cleanup();
  timer_done();
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

  for(i=CC_MIN; i <= CC_MAX; i++) psio_open(i,0);
}

void title(void)
{
  fprintf(outfile, "\n");
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t*         CCSORT         *\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\n");
}

void exit_io(void)
{
  int i;
  for(i=CC_MIN; i <= CC_MAX; i++)  psio_close(i,1);
  psio_done();
  ip_done();
  tstop(outfile);
  fclose(infile);
  fclose(outfile);
}

char *gprgid()
{
   char *prgid = "CCSORT";

   return(prgid);
}

void init_ioff(void)
{
  int i;
  ioff = init_int_array(IOFF_MAX);
  ioff[0] = 0;
  for(i=1; i < IOFF_MAX; i++) ioff[i] = ioff[i-1] + i;
}

void cleanup(void)
{
  int i;
  PSI_FPTR next;

  psio_write_entry(CC_INFO, "Reference Energy", (char *) &(moinfo.eref),
		   sizeof(double));
  
  free(moinfo.orbspi);
  free(moinfo.clsdpi);
  free(moinfo.openpi);
  free(moinfo.uoccpi);
  free(moinfo.fruocc);
  free(moinfo.frdocc);
  for(i=0; i < moinfo.nirreps; i++)
      free(moinfo.labels[i]);
  free(moinfo.labels);
  free(moinfo.occ);
  free(moinfo.vir);
  free(moinfo.socc);
  free(moinfo.all_occ);
  free(moinfo.all_vir);
  free(moinfo.all_socc);
  free(moinfo.frozen);
  free(moinfo.cc_occ);
  free(moinfo.cc_vir);
  free(moinfo.qt_occ);
  free(moinfo.qt_vir);
  free(moinfo.occ_sym);
  free(moinfo.vir_sym);
  free(moinfo.cc_allocc);
  free(moinfo.cc_allvir);
  free(moinfo.qt_allocc);
  free(moinfo.qt_allvir);
  free(moinfo.allocc_sym);
  free(moinfo.allvir_sym);
  free(moinfo.occ_off);
  free(moinfo.vir_off);
  free(moinfo.all_occ_off);
  free(moinfo.all_vir_off);
  free(moinfo.occpi);
  free(moinfo.virtpi);
  free(moinfo.all_occpi);
  free(moinfo.all_virtpi);

  free(ioff);

}
