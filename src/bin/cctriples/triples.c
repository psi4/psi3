#include <stdio.h>
#include <stdlib.h>
#include <ip_libv1.h>
#include <libciomr.h>
#include <psio.h>
#include <dpd.h>
#include <qt.h>
#include "globals.h"

void init_io(void);
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
void setup(void);
int **cacheprep(int level, int *cachefiles);

int main(int argc, char *argv[])
{
  int refwfn, errcod;
  char *junk;
  double ETAAA, ETAAB, ETABB, ETBBB, ET;
  long int memory;
  int **cachelist, *cachefiles;
  
  init_io();
  title();

  timer_init();
  timer_on("CCtriples");

  get_moinfo();
  fndcor(&(memory),infile,outfile);

  cachefiles = init_int_array(PSIO_MAXUNIT);
  cachelist = cacheprep(2, cachefiles);

  dpd_init(0, moinfo.nirreps, memory, 0, cachefiles, cachelist, NULL,
	   2, moinfo.occpi, moinfo.occ_sym, moinfo.virtpi, moinfo.vir_sym);

  errcod = ip_string("REFERENCE", &(junk), 0);
  if(!strcmp(junk, "RHF")) refwfn = 0;
  else if(!strcmp(junk, "ROHF")) refwfn = 1;
  else {
    printf("Invalid value of input keyword REFERENCE: %s\n", junk);
    exit(2);
  }

  if(refwfn == 0) {
    /*
      ET = ET_RHF_noddy();
      fprintf(outfile, "\t(T) energy                    = %20.15f\n", ET);
      fprintf(outfile, "\tTotal CCSD(T) energy          = %20.15f\n", 
      ET + moinfo.ecc + moinfo.eref);
    */

    /* New matrix-based code */
    ET = ET_RHF();
    fprintf(outfile, "\t(T) energy                    = %20.15f\n", ET);
    fprintf(outfile, "\tTotal CCSD(T) energy          = %20.15f\n", 
	    ET + moinfo.ecc + moinfo.eref);
    ET = ET_RHF_allijk();
    fprintf(outfile, "\t(T) energy                    = %20.15f\n", ET);
    fprintf(outfile, "\tTotal CCSD(T) energy          = %20.15f\n", 
	    ET + moinfo.ecc + moinfo.eref);
  }
  else if(refwfn == 1) {
    /*    setup(); */

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

  dpd_close(0);
  cleanup();

  timer_off("CCtriples");
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
  ip_done();
  tstop(outfile);
  fclose(infile);
  fclose(outfile);
}

char *gprgid()
{
   char *prgid = "CCTRIPLES";

   return(prgid);
}
