/*
**  CCDENSITY: Program to calculate the coupled-cluster one- and 
**             two-particle densities.
*/

#include <stdio.h>
#include <stdlib.h>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libiwl/iwl.h>
#include <math.h>
#include "globals.h"

/* Function prototypes */
void init_io(void);
void title(void);
void get_moinfo(void);
void get_frozen(void);
void get_params(void);
void exit_io(void);
void onepdm(void);
void sortone(void);
void twopdm(void);
void energy(void);
void resort_tei(void);
void resort_gamma(void);
void lag(void);
void build_X(void);
void build_A(void);
void build_Z(void);
void relax_I(void);
void relax_D(void);
void sortI(void);
void fold(void);
void deanti(void);
void add_ref(struct iwlbuf *OutBuf);
void add_core(struct iwlbuf *OutBuf);
void dump(struct iwlbuf *OutBuf);
void kinetic(void);
void probable(void);
int **cacheprep_rhf(int level, int *cachefiles);
int **cacheprep_uhf(int level, int *cachefiles);
void cachedone_rhf(int **cachelist);
void cachedone_uhf(int **cachelist);

int main(int argc, char *argv[])
{
  int **cachelist, *cachefiles;
  struct iwlbuf OutBuf;
  
  init_io();
  title();
  get_moinfo();
  get_frozen();
  get_params();

  if(moinfo.nfzc || moinfo.nfzv) {
    fprintf(outfile, "\n\tGradients involving frozen orbitals not yet available.\n");
    exit(2);
  }

  cachefiles = init_int_array(PSIO_MAXUNIT);

  if(params.ref == 0 || params.ref == 1) { /** RHF or ROHF **/
    cachelist = cacheprep_rhf(params.cachelev, cachefiles);

    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL, 
	     2, moinfo.occpi, moinfo.occ_sym, moinfo.virtpi, moinfo.vir_sym);

  }
  else if(params.ref == 2) { /** UHF **/
    cachelist = cacheprep_uhf(params.cachelev, cachefiles);

    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, 
	     cachelist, NULL, 4, moinfo.aoccpi, moinfo.aocc_sym, moinfo.avirtpi,
	     moinfo.avir_sym, moinfo.boccpi, moinfo.bocc_sym, moinfo.bvirtpi, moinfo.bvir_sym);
  }

  onepdm();
  twopdm();
  if(!params.aobasis) energy();
  sortone();
  kinetic();

  /*

  dpd_init(1, moinfo.nirreps, params.memory, 2, frozen.occpi, frozen.occ_sym,
  frozen.virtpi, frozen.vir_sym);
  */

  /*  if(moinfo.nfzc || moinfo.nfzv) {
      resort_gamma();
      resort_tei();
      } */

  lag();
  build_X();
  build_A();
  build_Z();

  /*  
      dpd_close(0);
      dpd_close(1);
  */


  relax_I();
  relax_D();
  sortone();
  sortI();
  fold();
  deanti();

  iwl_buf_init(&OutBuf, params.tpdmfile, params.tolerance, 0, 0);

  add_core(&OutBuf);
  add_ref(&OutBuf);
  dump(&OutBuf);
  iwl_buf_flush(&OutBuf, 1);
  iwl_buf_close(&OutBuf, 1);

  dpd_close(0);

  if(params.ref == 2) cachedone_uhf(cachelist);
  else cachedone_rhf(cachelist);
  free(cachefiles);

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
  fprintf(outfile, "\t\t\t*        CCDENSITY       *\n");
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
   char *prgid = "CCDENSITY";

   return(prgid);
}
