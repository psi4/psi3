/*
**  CCHBAR: Program to calculate the elements of the CCSD HBAR matrix.
*/

#include <stdio.h>
#include <stdlib.h>
#include <ip_libv1.h>
#include <psio.h>
#include <libciomr.h>
#include <dpd.h>
#include "globals.h"

/* Set this flag to 1 if you want re-run the code without re-running 
   ccenergy for testing purposes.  This will prevent CCLAMBDA from
   over-writing the incomplete WMBEJ HBAR matrix elements computed and
   used by CCENERGY.  CCLAMBDA will give incorrect results
   unless this flag is set to 0. */
#define REDO 0

/* Function prototypes */
void init_io(void);
void title(void);
void get_moinfo(void);
void get_params(void);
void exit_io(void);
void F_build(void);
void Wmbej_build(void);
void Wmnie_build(void);
void Wmbij_build(void);
void Wamef_build(void);
void Wabei_build(void);
void purge(void);
void cleanup(void);
int **cacheprep(int level, int *cachefiles);

int main(int argc, char *argv[])
{
  int **cachelist, *cachefiles;

  init_io();
  title();
  get_moinfo();
  get_params();

  cachefiles = init_int_array(PSIO_MAXUNIT);
  cachelist = cacheprep(params.cachelev, cachefiles);

  dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL,
           2, moinfo.occpi, moinfo.occ_sym, moinfo.virtpi, moinfo.vir_sym);

  if(!REDO) F_build();
  if(!REDO) Wmbej_build();
  Wamef_build();
  Wmnie_build();
  Wmbij_build();
  Wabei_build();

  purge();
  dpd_close(0);
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
  fprintf(outfile, "\t\t\t*         CCHBAR         *\n");
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
   char *prgid = "CCHBAR";

   return(prgid);
}
