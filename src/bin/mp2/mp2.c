#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <file30.h>
#include <qt.h>
#include <psio.h>
#include "MOInfo.h"
#include "Params.h"
#include "globals.h"

/* First definitions of globals */
FILE *infile, *outfile;
int *ioff;
struct MOInfo moinfo;
struct Params params;

/* Max length of ioff array */
#define IOFF 32641

/* Function Prototypes */
void init_io();
void title();
void init_ioff();
void get_parameters();
void get_moinfo();
void energy();
void exit_io();

int main()
{
  init_io();
  title();
  get_parameters();
  get_moinfo();
  init_ioff();
  energy();
  exit_io();
  exit(0);
}

void init_io(void)
{
  ffile(&infile,"input.dat",2);
  ffile(&outfile,"output.dat",1);
  tstart(outfile);
  ip_set_uppercase(1);
  ip_initialize(infile,outfile);
  ip_cwk_clear();
  ip_cwk_add(":DEFAULT");
  ip_cwk_add(":MP2");
  
  psio_init();
  file30_init();
}

void title(void)
{
  fprintf(outfile, "\n");
  fprintf(outfile,"\t***************************************************\n");
  fprintf(outfile,"\t*  MP2: Program to determine the MP2 Energy for   *\n");
  fprintf(outfile,"\t*               closed-shell systems.             *\n");
  fprintf(outfile,"\t*                                                 *\n");
  fprintf(outfile,"\t*                Daniel Crawford                  *\n");
  fprintf(outfile,"\t*                September  1995                  *\n");
  fprintf(outfile,"\t***************************************************\n");
  fprintf(outfile, "\n");
}

void init_ioff(void)
{
  int i;
  ioff = (int *) malloc(IOFF * sizeof(int));
  if(ioff == NULL) {
      fprintf(stderr, "(mp2): error malloc'ing ioff array\n");
      exit(0);
          }
  ioff[0] = 0;
  for(i=1; i < IOFF; i++) {
      ioff[i] = ioff[i-1] + i;
    }
}

void exit_io(void)
{
  tstop(outfile);
  ip_done();
  psio_done();
  file30_close();
}

void get_parameters(void)
{
  int errcod, tol;

  params.print_lvl = 0;
  errcod = ip_data("PRINT_LVL", "%d", &(params.print_lvl),0);
  errcod = ip_string("WFN", &(params.wfn),0);
  if(errcod == IPE_KEY_NOT_FOUND) {
      params.wfn = (char *) malloc(sizeof(char)*5);
      strcpy(params.wfn, "MP2");
    }
  params.nfile = 72;
  errcod = ip_data("NFILE", "%d", &(params.nfile),0);
  params.tolerance = 1e-14;
  errcod = ip_data("TOLERANCE","%d",&(tol),0);
  if(errcod == IPE_OK) {
      params.tolerance = 1.0*pow(10.0,(double) -tol);
    }
  params.keep_integrals=0;
  errcod = ip_boolean("KEEP_INTEGRALS", &(params.keep_integrals),0);

  fprintf(outfile,"\tInput Parameters:\n");
  fprintf(outfile,"\t-----------------\n");
  fprintf(outfile,"\tWavefunction           =  %s\n", params.wfn);
  fprintf(outfile,"\tPrint Level            =  %d\n", params.print_lvl);
  fprintf(outfile,"\tN File                 =  %d\n", params.nfile);
  fprintf(outfile,"\tTolerance              =  %3.1e\n", params.tolerance);
  fprintf(outfile,"\tKeep Integrals         =  %s\n", (params.keep_integrals ? "Yes": "No"));
}

void get_moinfo(void)
{
  int i,errcod,size;

  moinfo.frdocc = init_int_array(8);
  errcod = ip_count("FROZEN_DOCC",&size,0);
  if(errcod == IPE_OK) {
      for(i=0; i < size; i++)
          errcod = ip_data("FROZEN_DOCC","%d",(&(moinfo.frdocc)[i]),1,i);
    }
  moinfo.fruocc = init_int_array(8);
  errcod = ip_count("FROZEN_UOCC",&size,0);
  if(errcod == IPE_OK) {
      for(i=0; i < size; i++)
          errcod = ip_data("FROZEN_UOCC","%d",(&(moinfo.fruocc)[i]),1,i);
    }
  
  moinfo.nmo = file30_rd_nmo();
  moinfo.nirreps = file30_rd_nirreps();
  moinfo.iopen = file30_rd_iopen();
  moinfo.labels = file30_rd_irr_labs();
  moinfo.orbspi = file30_rd_orbspi();
  moinfo.clsdpi = file30_rd_clsdpi();
  moinfo.openpi = file30_rd_openpi();
  moinfo.enuc = file30_rd_enuc();
  moinfo.escf = file30_rd_escf();
  moinfo.evals = file30_rd_evals();

  
  moinfo.virtpi = init_int_array(moinfo.nirreps);
  for(i=0; i < moinfo.nirreps; i++) {
      moinfo.virtpi[i] = moinfo.orbspi[i]-moinfo.clsdpi[i]-moinfo.openpi[i];
    }

  fprintf(outfile,"\n\tFile30 Parameters:\n");
  fprintf(outfile,"\t------------------\n");
  fprintf(outfile,"\tNumber of irreps = %d\n",moinfo.nirreps);
  fprintf(outfile,"\tNumber of MOs    = %d\n\n",moinfo.nmo);
  fprintf(outfile,"\tLabel\t# MOs\t# FZDC\t# DOCC\t# SOCC\t# VIRT\t# FZVR\n");
  fprintf(outfile,"\t-----\t-----\t------\t------\t------\t------\t------\n");
  for(i=0; i < moinfo.nirreps; i++) {
      fprintf(outfile,
              "\t %s\t   %d\t    %d\t    %d\t    %d\t    %d\t    %d\n",
              moinfo.labels[i],moinfo.orbspi[i],moinfo.frdocc[i],
              moinfo.clsdpi[i],moinfo.openpi[i],moinfo.virtpi[i],
              moinfo.fruocc[i]);
    }
  fprintf(outfile,"\n\tNuclear Repulsion Energy    = %20.10f\n",moinfo.enuc);
  fprintf(outfile,  "\tTotal SCF Energy            = %20.10f\n",moinfo.escf);

  /* Set up some other useful parameters */
  moinfo.noeints = moinfo.nmo*(moinfo.nmo+1)/2;
  moinfo.nteints = moinfo.noeints*(moinfo.noeints+1)/2;
}
