#include <stdio.h>
#include <stdlib.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include "globals.h"

void init_io(int argc, char *argv[]);
void title(void);
void get_moinfo(void);
void exit_io(void);
void cleanup(void);
double ET_RHF(void);
double ET_RHF_allijk(void);
double ET_RHF_noddy(void);
double ET_AAA(void);
double ET_AAB(void);
double ET_ABB(void);
double ET_BBB(void);
double ET_UHF_AAA(void);
double ET_UHF_BBB(void);
double ET_UHF_AAB(void);
double ET_UHF_ABB(void);
double ET_UHF_AAA_noddy(void);
double ET_UHF_BBB_noddy(void);
double ET_UHF_AAB_noddy(void);
double ET_UHF_ABB_noddy(void);
void count_ijk(void);
void setup(void);
int **cacheprep_rhf(int level, int *cachefiles);
int **cacheprep_uhf(int level, int *cachefiles);
void cachedone_uhf(int **cachelist);
void cachedone_rhf(int **cachelist);

int main(int argc, char *argv[])
{
  double ETAAA, ETAAB, ETABB, ETBBB, ET;
  long int memory;
  int **cachelist, *cachefiles;
  dpdfile2 T1;
  double **geom, *zvals, value;
  FILE *efile;
  int i, errcod, natom;
  
  init_io(argc, argv);
  title();

  timer_init();
  timer_on("CCtriples");

  get_moinfo();
  fndcor(&(memory),infile,outfile);

  errcod = ip_string("WFN", &(params.wfn), 0);
  if(strcmp(params.wfn, "CCSD") && strcmp(params.wfn, "CCSD_T")) {
    fprintf(outfile, "Invalid value of input keyword WFN: %s\n", params.wfn);
    exit(2);
  }

  cachefiles = init_int_array(PSIO_MAXUNIT);

  if(params.ref == 0) { /*** RHF ***/
    cachelist = cacheprep_rhf(2, cachefiles);

    dpd_init(0, moinfo.nirreps, memory, 0, cachefiles, cachelist, NULL,
	     2, moinfo.occpi, moinfo.occ_sym, moinfo.virtpi, moinfo.vir_sym);
  }
  else if(params.ref == 2) { /*** UHF ***/
    cachelist = cacheprep_uhf(2, cachefiles);

    dpd_init(0, moinfo.nirreps, memory, 0, cachefiles, 
	     cachelist, NULL, 4, moinfo.aoccpi, moinfo.aocc_sym, moinfo.avirtpi,
	     moinfo.avir_sym, moinfo.boccpi, moinfo.bocc_sym, moinfo.bvirtpi, moinfo.bvir_sym);
  }

  count_ijk();
  fflush(outfile);

  if(params.ref == 0) { /** RHF **/

    ET = ET_RHF();
    fprintf(outfile, "\t(T) energy                    = %20.15f\n", ET);
    fprintf(outfile, "\tTotal CCSD(T) energy          = %20.15f\n", 
	    ET + moinfo.ecc + moinfo.eref);
  }
  else if(params.ref == 1) { /** ROHF --- don't use this right now! **/

    fprintf(outfile, "\nROHF-CCSD(T) is not yet available...\n");
    exit(2); 

    ETAAA = ET_AAA();
    fprintf(outfile, "\tAAA (T) energy                = %20.15f\n", ETAAA);
    ETAAB = ET_AAB();
    fprintf(outfile, "\tAAB (T) energy                = %20.15f\n", ETAAB);
    ETABB = ET_ABB();
    fprintf(outfile, "\tABB (T) energy                = %20.15f\n", ETABB);
    ETBBB = ET_BBB();
    fprintf(outfile, "\tBBB (T) energy                = %20.15f\n", ETBBB);
    ET = ETAAA + ETAAB + ETABB + ETBBB;
    fprintf(outfile, "\t(T) energy                    = %20.15f\n", ET);
    fprintf(outfile, "\tTotal CCSD(T) energy          = %20.15f\n", 
	    ET + moinfo.ecc + moinfo.eref);
  }
  else if(params.ref == 2) { /** UHF **/

    ETAAA = ET_UHF_AAA();
    fprintf(outfile, "\tAAA (T) energy                = %20.15f\n", ETAAA);
    fflush(outfile);

    ETBBB = ET_UHF_BBB();
    fprintf(outfile, "\tBBB (T) energy                = %20.15f\n", ETBBB);
    fflush(outfile);

    ETAAB = ET_UHF_AAB();
    fprintf(outfile, "\tAAB (T) energy                = %20.15f\n", ETAAB);
    fflush(outfile);

    ETABB = ET_UHF_ABB();
    fprintf(outfile, "\tABB (T) energy                = %20.15f\n", ETABB);
    fflush(outfile);

    ET = ETAAA + ETAAB + ETABB + ETBBB;
    fprintf(outfile, "\t(T) energy                    = %20.15f\n", ET);
    fprintf(outfile, "\tTotal CCSD(T) energy          = %20.15f\n", 
	    ET + moinfo.ecc + moinfo.eref);
  }

  fprintf(outfile, "\n");

  /* Write total energy to the checkpoint file */
  chkpt_init(PSIO_OPEN_OLD);
  chkpt_wt_etot(ET+moinfo.ecc+moinfo.eref);
  chkpt_close();

  /* Write pertinent data to energy.dat */
  if(!strcmp(params.wfn,"CCSD_T")) {
    chkpt_init(PSIO_OPEN_OLD);
    natom = chkpt_rd_natom();
    geom = chkpt_rd_geom();
    zvals = chkpt_rd_zvals();
    chkpt_close();
    ffile(&efile,"energy.dat",1);
    fprintf(efile, "*\n");
    for(i=0; i < natom; i++) 
      fprintf(efile, " %4d   %5.2f     %13.10f    %13.10f    %13.10f\n",
	      i+1, zvals[i], geom[i][0], geom[i][1], geom[i][2]);
    free_block(geom);  free(zvals);
    fprintf(efile, "SCF(30)   %22.12f\n", moinfo.escf);
    fprintf(efile, "REF(100)  %22.12f\n", moinfo.eref);
    fprintf(efile, "CCSD      %22.12f\n", (moinfo.ecc+moinfo.eref));
    fprintf(efile, "CCSD(T)   %22.12f\n", (ET+ moinfo.ecc+moinfo.eref));
    fclose(efile);
  }

  dpd_close(0);

  if(params.ref == 2) cachedone_uhf(cachelist);
  else cachedone_rhf(cachelist);
  free(cachefiles);

  cleanup();

  timer_off("CCtriples");
  timer_done();

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

  psi_start(argc-1,argv+1,0);
  ip_cwk_add(progid);
  free(progid);
  tstart(outfile);

  psio_init();
  for(i=CC_MIN; i <= CC_MAX; i++) psio_open(i,1);
}

void title(void)
{
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t*        CCTRIPLES       *\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t**************************\n");
}

void exit_io(void)
{
  int i;
  for(i=CC_MIN; i <= CC_MAX; i++) psio_close(i,1);
  psio_done();
  tstop(outfile);
  psi_stop();
}

char *gprgid()
{
   char *prgid = "CCTRIPLES";

   return(prgid);
}
