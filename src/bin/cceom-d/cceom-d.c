/*
**  CCEOM-D : Compute the non-R0 EOM parts of the CC 1pdm and 2pdm
*/

#include <stdio.h>
#include <stdlib.h>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libiwl/iwl.h>
#include <math.h>
#include <psifiles.h>
#include "globals.h"

/* Function prototypes */
void init_io(int argc, char *argv[]);
void title(void);
void get_moinfo(void);
void get_params(void);
void exit_io(void);
void onepdm(void);
void energy(void);
int **cacheprep_rhf(int level, int *cachefiles);
int **cacheprep_uhf(int level, int *cachefiles);
void cachedone_rhf(int **cachelist);
void cachedone_uhf(int **cachelist);
void norm_Dpq(void);
void intermediates(int R_irr);

int main(int argc, char *argv[])
{
  int **cachelist, *cachefiles, R_irr = 0;
  struct iwlbuf OutBuf;
  struct iwlbuf OutBuf_AA, OutBuf_BB, OutBuf_AB;
  
  init_io(argc,argv);
  title();
  get_moinfo();
  get_params();

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

  intermediates(R_irr);
  onepdm();
  norm_Dpq(); /* for debugging */
  /* twopdm(); */
  if(!params.aobasis) energy();

  dpd_close(0);

  if(params.ref == 2) cachedone_uhf(cachelist);
  else cachedone_rhf(cachelist);
  free(cachefiles);

  cleanup(); 
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

  psi_start(argc-1,argv+1,0);
  ip_cwk_add(progid);
  free(progid);
  tstart(outfile);

  psio_init();

  /* Open all dpd data files here */
  for(i=CC_MIN; i <= CC_MAX; i++) psio_open(i,1);
}

void title(void)
{
  fprintf(outfile, "\n");
  fprintf(outfile, "\t\t\t*****************************\n");
  fprintf(outfile, "\t\t\t*          CCEOM-D          *\n");
  fprintf(outfile, "\t\t\t*  By Rollin King    2003   *\n");
  fprintf(outfile, "\t\t\t*****************************\n");
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
   char *prgid = "CCDENSITY";

   return(prgid);
}

void norm_Dpq(void) {
  double nDIA, nDia, nDIJ, nDij, nDAI, nDai, nDAB, nDab;
  dpdfile2 DIA, Dia, DIJ, Dij, DAI, Dai, DAB, Dab;

  dpd_file2_init(&DIJ, EOM_D, 0, 0, 0, "DIJ");
  nDIJ = dpd_file2_dot_self(&DIJ);
  dpd_file2_close(&DIJ);

  dpd_file2_init(&Dij, EOM_D, 0, 0, 0, "Dij");
  nDij = dpd_file2_dot_self(&Dij);
  dpd_file2_close(&Dij);

  dpd_file2_init(&DAB, EOM_D, 0, 1, 1, "DAB");
  nDAB = dpd_file2_dot_self(&DAB);
  dpd_file2_close(&DAB);

  dpd_file2_init(&Dab, EOM_D, 0, 1, 1, "Dab");
  nDab = dpd_file2_dot_self(&Dab);
  dpd_file2_close(&Dab);

  dpd_file2_init(&DAI, EOM_D, 0, 0, 1, "DAI");
  nDAI = dpd_file2_dot_self(&DAI);
  dpd_file2_close(&DAI);

  dpd_file2_init(&Dai, EOM_D, 0, 0, 1, "Dai");
  nDai = dpd_file2_dot_self(&Dai);
  dpd_file2_close(&Dai);

  dpd_file2_init(&DIA, EOM_D, 0, 0, 1, "DIA");
  nDIA = dpd_file2_dot_self(&DIA);
  dpd_file2_close(&DIA);

  dpd_file2_init(&Dia, EOM_D, 0, 0, 1, "Dia");
  nDia = dpd_file2_dot_self(&Dia);
  dpd_file2_close(&Dia);

  fprintf(outfile,"Self-dots of Dpq * Dpq \n");

  fprintf(outfile,"<DIJ|DIJ> = %15.10lf\n", nDIJ);
  fprintf(outfile,"<Dij|Dij> = %15.10lf\n", nDij);

  fprintf(outfile,"<DAB|DAB> = %15.10lf\n", nDAB);
  fprintf(outfile,"<Dab|Dab> = %15.10lf\n", nDab);

  fprintf(outfile,"<DAI|DAI> = %15.10lf\n", nDAI);
  fprintf(outfile,"<Dai|Dai> = %15.10lf\n", nDai);

  fprintf(outfile,"<DIA|DIA> = %15.10lf\n", nDIA);
  fprintf(outfile,"<Dia|Dia> = %15.10lf\n", nDia);

  return;
}
