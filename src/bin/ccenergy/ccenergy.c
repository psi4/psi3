/*
**  CCENERGY: Program to calculate coupled cluster energies.
*/

#include <stdio.h>
#include <stdlib.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libfile30/file30.h>
#include <libqt/qt.h>
#include "globals.h"

/* Function prototypes */
void init_io(void);
void title(void);
void get_moinfo(void);
void get_params(void);
void init_amps(void);
void tau_build(void);
void taut_build(void);
double energy(void);
void sort_amps(void);
void Fae_build(void);
void Fmi_build(void);
void Fme_build(void);
void t1_build(void);
void Wmnij_build(void);
void Z_build(void);
void Y_build(void);
void X_build(void);
void Wmbej_build(void);
void t2_build(void);
void tsave(void);
int converged(void);
double diagnostic(void);
double d1diag(void);
void exit_io(void);
void cleanup(void);
void update(void);
void diis(int iter);
void ccdump(void);
int **cacheprep_uhf(int level, int *cachefiles);
int **cacheprep_rhf(int level, int *cachefiles);
void cachedone_rhf(int **cachelist);
void cachedone_uhf(int **cachelist);
void memchk(void);
struct dpd_file4_cache_entry *priority_list(void);
void spinad_amps(void);
void status(char *, FILE *);

/* local correlation functions */
void local_init(void);
void local_done(void);

int main(int argc, char *argv[])
{
  int done=0;
  int h, i, j, a, b, row, col, natom;
  double **geom, *zvals, value;
  FILE *efile;
  int **cachelist, *cachefiles;
  struct dpd_file4_cache_entry *priority;
  dpdfile2 T1;

  moinfo.iter=0;
  
  init_io();
  title();

  timer_init();
  timer_on("CCEnergy");
  
  get_moinfo();
  get_params();

  cachefiles = init_int_array(PSIO_MAXUNIT);

  if(params.ref == 2) { /** UHF **/
    cachelist = cacheprep_uhf(params.cachelev, cachefiles);

    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, 
	     cachelist, NULL, 4, moinfo.aoccpi, moinfo.aocc_sym, moinfo.avirtpi,
	     moinfo.avir_sym, moinfo.boccpi, moinfo.bocc_sym, moinfo.bvirtpi, moinfo.bvir_sym);

    if(params.aobasis) { /* Set up new DPD's for AO-basis algorithm */
      dpd_init(1, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL, 
               4, moinfo.aoccpi, moinfo.aocc_sym, moinfo.orbspi, moinfo.orbsym, 
	       moinfo.boccpi, moinfo.bocc_sym, moinfo.orbspi, moinfo.orbsym);
      dpd_set_default(0);
    }

  }
  else { /** RHF or ROHF **/
    cachelist = cacheprep_rhf(params.cachelev, cachefiles);

    priority = priority_list();

    dpd_init(0, moinfo.nirreps, params.memory, params.cachetype, cachefiles, 
	     cachelist, priority, 2, moinfo.occpi, moinfo.occ_sym, 
	     moinfo.virtpi, moinfo.vir_sym);
   
    if(params.aobasis) { /* Set up new DPD for AO-basis algorithm */
      dpd_init(1, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL, 
	       2, moinfo.occpi, moinfo.occ_sym, moinfo.orbspi, moinfo.orbsym);
      dpd_set_default(0);
    }

  }

  if(params.local) local_init(); 
 
  init_amps();
  tau_build();
  taut_build();
  fprintf(outfile, "\t                     Solving CCSD Equations\n");
  fprintf(outfile, "\t                     ----------------------\n");
  fprintf(outfile, "\tIter             Energy               RMS       T1Diag      D1Diag\n");
  fprintf(outfile, "\t----     ---------------------     --------   ----------  ----------\n");
  moinfo.ecc = energy();
  moinfo.t1diag = diagnostic();
  moinfo.d1diag = d1diag();
  update();
  for(moinfo.iter=1; moinfo.iter <= params.maxiter; moinfo.iter++) {

    timer_on("sort_amps");
    sort_amps();
    timer_off("sort_amps");

    timer_on("F build");
    Fme_build(); Fae_build(); Fmi_build();
    if(params.print & 2) status("F intermediates", outfile);
    timer_off("F build");

    timer_on("T1 Build");
    t1_build();
    if(params.print & 2) status("T1 amplitudes", outfile);
    timer_off("T1 Build");

    timer_on("Wmbej build");
    Wmbej_build();
    if(params.print & 2) status("Wmbej", outfile);
    timer_off("Wmbej build");

    Z_build();
    if(params.print & 2) status("Z", outfile);
    Wmnij_build();
    if(params.print & 2) status("Wmnij", outfile);

    timer_on("T2 Build");
    t2_build();
    if(params.print & 2) status("T2 amplitudes", outfile);
    timer_off("T2 Build");

    if(converged()) {
      done = 1;
      tsave();
      tau_build(); taut_build();
      moinfo.ecc = energy();
      moinfo.t1diag = diagnostic();
      moinfo.d1diag = d1diag();
      sort_amps();
      update();
      fprintf(outfile, "\n\tIterations converged.\n");
      fflush(outfile);
      break;
    }
    if(params.diis) diis(moinfo.iter);
    tsave();
    tau_build(); taut_build();
    moinfo.ecc = energy();
    moinfo.t1diag = diagnostic();
    moinfo.d1diag = d1diag();
    update();
  }
  fprintf(outfile, "\n");
  if(!done) {
    fprintf(outfile, "\t ** Wave function not converged to %2.1e ** \n",
	    params.convergence);
    fflush(outfile);
    if(params.aobasis) dpd_close(1);
    dpd_close(0);
    cleanup();
    timer_off("CCEnergy");
    timer_done();
    exit_io();
    exit(1);
  }

  fprintf(outfile, "\tSCF energy       (file30)  = %20.15f\n", moinfo.escf);
  fprintf(outfile, "\tReference energy (file100) = %20.15f\n", moinfo.eref);
  fprintf(outfile, "\tCCSD correlation energy    = %20.15f\n", moinfo.ecc);
  fprintf(outfile, "\tTotal CCSD energy          = %20.15f\n", 
          moinfo.eref + moinfo.ecc);
  fprintf(outfile, "\n");

  /* Write pertinent data to energy.dat for Dr. Yamaguchi */
  if(!strcmp(params.wfn,"CCSD")) {
    file30_init();
    natom = file30_rd_natom();
    geom = file30_rd_geom();
    zvals = file30_rd_zvals();
    file30_close();
    ffile(&efile, "energy.dat",1);
    fprintf(efile, "*\n");
    for(i=0; i < natom; i++) 
      fprintf(efile, " %4d   %5.2f     %13.10f    %13.10f    %13.10f\n",
	      i+1, zvals[i], geom[i][0], geom[i][1], geom[i][2]);
    free_block(geom);  free(zvals);
    fprintf(efile, "SCF(30)   %22.12f\n", moinfo.escf);
    fprintf(efile, "REF(100)  %22.12f\n", moinfo.eref);
    fprintf(efile, "CCSD      %22.12f\n", (moinfo.ecc+moinfo.eref));
    fclose(efile);
  }

  /* Generate the spin-adapted RHF amplitudes for later codes */
  if(params.ref == 0) {
    timer_on("spinad Amps");
    spinad_amps();
    timer_off("spinad Amps");
  }

  if(params.aobasis) dpd_close(1);
  dpd_close(0);

  if(params.ref == 2) cachedone_uhf(cachelist);
  else cachedone_rhf(cachelist);
  free(cachefiles);

  if(params.local) local_done();
  
  cleanup();

  timer_off("CCEnergy");
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
  for(i=CC_MIN; i <= CC_LAMPS; i++) psio_open(i,1);
  for(i=CC_HBAR; i <= CC_MAX; i++) psio_open(i,0);
}

void title(void)
{
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t*        CCENERGY        *\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t**************************\n");
}

void exit_io(void)
{
  int i;
  for(i=CC_MIN; i <= CC_MAX; i++) psio_close(i,1);
  psio_done();

  free_ptrs();
  ip_done();
  tstop(outfile);
  fclose(infile);
  fclose(outfile);
}

char *gprgid()
{
   char *prgid = "CCENERGY";

   return(prgid);
}

void memchk(void)
{
  pid_t mypid;
  FILE *memdat;
  char comm[80];

  mypid = getpid();

  memdat = freopen("output.dat","a",stdout);
 
  sprintf(comm, "grep \"VmSize\" /proc/%d/status", (int) mypid);

  system(comm);
}
